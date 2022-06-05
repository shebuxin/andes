from andes.interop.pandapower import to_pandapower
from andes.interop.pandapower import make_GSF, build_group_table
import gurobipy as gb
import pandas as pd
import numpy as np
import logging
import torch
import os
logger = logging.getLogger(__name__)

from opf import dcopf # base class

'''
    Vittual Inertia Scheduling (vis) solves economic dispatch problem
    consiering the dynamic frequency constriants and vsg invertia 
    support reserve.

    The base class is inherited from 'dcopf' in opf.py

    1) vis1: dcopf + fnadir/RoCof (ML linearization)

    2) vis2: dcopf + fnadir/RoCof (ML linearization)
                   + VSG power reserve (ML linearization)

    3) vis3: dcopf + fnadir/RoCof (ML linearization)
                   + VSG power reserve (Final value theorem)


    Auxiliary functions: 

    1) loadnn: load network weights, bias and normalization parameters

'''

# ----------------------------- class: vis1 --------------------------------------
class vis1(dcopf):
    """
    vis1: dcopf + fnadir/RoCof (ML linearization)
    """

    def __init__(self, name='vis1', norm=None, nn=None, nn_num=64, dpe=0.01, rocof_lim = 0.01, nadir_lim=0.01):
        """
        Initialize high level parameters

        Input
        ---------
        name: str
            name
        norm: dict
            normalizatin parameter for fnadir and Pvsg prediction
        nn: dict
            neural network weight and bias
        nn_num: integer
            number of MLP nuerols (assume single layer MLP)
        dpe: float
            delta Pe, power mismatch, or load change
        """
        super().__init__(name = name)
        self.norm = norm 
        self.nn = nn
        self.nn_num = nn_num
        self.fnadir = nadir_lim
        self.rocof = rocof_lim
        self.dpe = dpe

    def from_andes(self, ssa, typeII=None, Sbase=100):
        """
        Initialize parameters from andes mdoel

        Input
        ---------
        ssa: andes model
        typeII: list
            idx of typeII generator (vsg inverter), i.g., ['PV_6', 'PV_7']
        Sbase: float
            system base power (MW)
        """

        super().from_andes(ssa)

        # Define typeII generatiors, defalt typeI (type=1)
        self.gen['type'] = 1
        if typeII:
            for idx in typeII:
                row = self.gen[self.gen['idx'] == idx].index[0]
                self.gen['type'].iloc[row] = 2

        # add new parameter
        self.scale = Sbase / self.mva
        self.gen['p_pre'] = 0
        self.gen['band'] = self.gen['pmax'] - self.gen['pmin']
        self.gen['Sn'] /= self.mva  # normalize Sn
        self.gen['K'] = 1

        # load new parameter from andes
        # Note: control parameters are normalized accroding to andes system base: Sbase = 100 MVA
        genrow = ssa.GENROU.as_df()
        regc = ssa.REGCV1.as_df()
        tgov = ssa.TGOV1N.as_df()

        tgov.rename(columns={'idx':'gov', 'syn':'idx'}, inplace=True)
        regc.rename(columns={'idx':'vsg', 'gen':'idx', 'M':'Mvsg', 'D': 'Dvsg'},  inplace=True)
        # merge dataframe
        genrow = pd.merge(left=genrow, right=tgov[['idx', 'R']], on='idx', how='left')
        genrow.rename(columns={'idx': 'syn', 'gen': 'idx'}, inplace=True)

        self.gen = pd.merge(left=self.gen, right=genrow[['idx', 'M','D', 'R']], on='idx', how='left')
        self.gen = pd.merge(left=self.gen, right=regc[['idx', 'Mvsg', 'Dvsg']], on='idx', how='left')
        # fill nan caused by merge
        self.gen.fillna(0, inplace=True)

        # norm control parameter according to Sbase
        self.gen['M'] /= self.scale
        self.gen['D'] /= self.scale
        self.gen['R'] *= self.scale
        self.gen['Mvsg'] /= self.scale
        self.gen['Dvsg'] /= self.scale

        self.update_dict()

    def build(self):
        # self.data_check()

        # --- build gurobi model ---
        self.update_dict()
        self.mdl = gb.Model(self.name)
        self._build_vars()
        self._build_obj()
        self._build_cons()

        logger.info('Successfully build vis0 model.')

    def _build_vars(self):
        GEN = self.gendict.keys()

        # --- uncontrollable generators limit to p0 ---
        gencp = self.gen.copy()
        gencp['pmax'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        gencp['pmin'][gencp.ctrl == 0] = gencp['p0'][gencp.ctrl == 0]
        # --- offline geenrators limit to 0 ---
        gencp['pmax'][gencp.u == 0] = 0
        gencp['pmin'][gencp.u == 0] = 0

        # --- gen: pg ---
        self.pg = self.mdl.addVars(GEN, name='pg', vtype=gb.GRB.CONTINUOUS, obj=0,
                              ub=gencp.pmax.tolist(), lb=gencp.pmin.tolist())
        # --- RegUp, RegDn --- !!! modify inverter reserve up and down in andes file, prd=0
        self.pru = self.mdl.addVars(GEN, name='pru', vtype=gb.GRB.CONTINUOUS, obj=0,
                               ub=gencp.band.tolist(), lb=[0] * gencp.shape[0])
        self.prd = self.mdl.addVars(GEN, name='prd', vtype=gb.GRB.CONTINUOUS, obj=0,
                               ub=gencp.band.tolist(), lb=[0] * gencp.shape[0])
        
        # --- Mvsg, Dvsg ---
        _, vsg = self._get_GENI_GENII_key()

        self.Mvsg = self.mdl.addVars(vsg, name='Mvsg', vtype=gb.GRB.CONTINUOUS, obj=0,
                               ub=[5]*len(vsg), lb=[0]*len(vsg))
        self.Dvsg = self.mdl.addVars(vsg, name='Dvsg', vtype=gb.GRB.CONTINUOUS, obj=0,
                               ub=[3]*len(vsg), lb=[0]*len(vsg))

        print('Successfully build var.')

    def _build_obj(self):
        GEN = self.gendict.keys()
        gendict = self.gendict
        costdict = self.costdict

        # --- minimize generation cost ---
        cost_pg = sum(
                        (self.pg[gen] * costdict[gen]['c1']
                        + self.pg[gen] * self.pg[gen] * costdict[gen]['c2']
                        + costdict[gen]['c0']) * gendict[gen]['u']  # online status
                        for gen in GEN
                    )

        # --- RegUp, RegDn cost ---
        cost_ru = sum(self.pru[gen] * costdict[gen]['cru'] * gendict[gen]['u'] for gen in GEN)
        cost_rd = sum(self.prd[gen] * costdict[gen]['crd'] * gendict[gen]['u'] for gen in GEN)
        cost_vsg = cost_ru + cost_rd

        self.obj = self.mdl.setObjective(expr=cost_pg + cost_vsg, sense=gb.GRB.MINIMIZE)
        print('Successfully build obj.')

    def _build_cons(self):
        # --- var idx ---
        ptotal = self.load.p0.sum()

        gendict = self.gendict
        linedict = self.linedict
        gen_gsfdict = self.gen_gsfdict

        GEN = gendict.keys()
        LINE = linedict.keys()

        # --- 01 power balance ---
        p_sum = sum(self.pg[gen] for gen in GEN)
        self.mdl.addConstr(p_sum == ptotal, name='PowerBalance')

        # --- 02 line limits ---
        for line in LINE:
            lhs1 = sum(self.pg[gen] * gen_gsfdict[gen][line] for gen in GEN)
            self.mdl.addConstr(lhs1+linedict[line]['sup'] <= linedict[line]['rate_a'], name=f'{line}_U')
            self.mdl.addConstr(lhs1+linedict[line]['sup'] >= -linedict[line]['rate_a'], name=f'{line}_D')

        # --- 03 dynamic frequency constraints --
        self._add_fcons()

        print('Successfully build cons.')

    def _add_fcons(self):

        gendict = self.gendict
        GENI, GENII = self._get_GENI_GENII_key()
        
        # --- Synthetic M/D/F/R ---
        Msys = sum(gendict[gen]['Sn'] * gendict[gen]['M'] for gen in GENI)
        Msys += sum(gendict[gen]['Sn'] * self.Mvsg[gen] for gen in GENII)
        Msys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Dsys = sum(gendict[gen]['Sn'] * gendict[gen]['D'] for gen in GENI)
        Dsys += sum(gendict[gen]['Sn'] * self.Dvsg[gen] for gen in GENII)
        Dsys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Rsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * gendict[gen]['Sn'] for gen in GENI)
        Rsys /= sum(gendict[gen]['Sn'] for gen in GENI)

        Fsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * self.pg[gen] for gen in GENI)
        Fsys /= sum(gendict[gen]['Sn'] for gen in GENI)
        
        # --- add gurobi var for fnadir ---
        self.af = self.mdl.addVars(self.nn_num, name='af', vtype=gb.GRB.BINARY)
        self.zf = self.mdl.addVars(self.nn_num, name='zf', vtype=gb.GRB.CONTINUOUS, obj=0, lb=[0]*self.nn_num)

        # --- add constraints ---
        # --- RoCof con----
        self.mdl.addConstr(self.rocof * Msys >= abs(self.dpe), name='RoCof')

        # --- fnadir con----
        # - self.fnadir <= fnadir <= self.fnadir
        fnorm = self.norm['fnorm']
        Msys_norm = (Msys - fnorm['M'].iloc[0]) / fnorm['M'].iloc[1] # iloc[0] is mean, iloc[1] is std
        Dsys_norm = (Dsys - fnorm['D'].iloc[0]) / fnorm['D'].iloc[1]
        Fsys_norm = (Fsys - fnorm['Fg'].iloc[0]) / fnorm['Fg'].iloc[1]
        Rsys_norm = (Rsys - fnorm['Rg'].iloc[0]) / fnorm['Rg'].iloc[1]

        hdown = -100
        hup = 100

        # nn hidden layer constraints
        zf_bar = []
        for i in range(self.nn_num):
            zf_bar_temp = Msys_norm * self.nn['fw1'][0].iloc[i] + \
                          Dsys_norm * self.nn['fw1'][1].iloc[i] + \
                          Fsys_norm * self.nn['fw1'][2].iloc[i] + \
                          Rsys_norm * self.nn['fw1'][3].iloc[i] + self.nn['fb1'][0].iloc[i]
            zf_bar.append(zf_bar_temp)

        self.mdl.addConstrs(self.zf[i] <= zf_bar[i] - hdown*(1-self.af[i]) for i in range(self.nn_num))
        self.mdl.addConstrs(self.zf[i] >= zf_bar[i] for i in range(self.nn_num))
        self.mdl.addConstrs(self.zf[i] <= hup * self.af[i] for i in range(self.nn_num))

        # nn output (fnadir) constraints
        fnadir_norm = sum(self.zf[i]*self.nn['fw2'][i].iloc[0] for i in range(self.nn_num)) + self.nn['fb2'].iloc[0]

        # denormalize fnadir
        fnadir_pred = fnadir_norm * fnorm['fnadir'].iloc[1] + fnorm['fnadir'].iloc[0]
        
        fnadir_pred *=self.dpe

        self.mdl.addConstr( fnadir_pred >= - self.fnadir, name=f'fnadir_D')
        self.mdl.addConstr( fnadir_pred <= self.fnadir, name=f'fnadir_U')


    def _get_GENI_GENII_key(self):
        gendict = self.gendict

        # --- filter Type II gen ---
        gendict_I, gendict_II= dict(), dict()
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 1:
                gendict_I[new_key] = new_value
        for (new_key, new_value) in gendict.items():
            if new_value['type'] == 2:
                gendict_II[new_key] = new_value

        return gendict_I.keys(), gendict_II.keys()

    def get_res(self):
        """
        Get resutlts, can be used after mdl.optimize().

        Returns
        -------
        DataFrame
            The output DataFrame contains setpoints ``pg``

        """
        self.build()
        self.mdl.optimize()

        status = self.mdl.status
        
        if status == gb.GRB.UNBOUNDED:
            print('The model cannot be solved because it is unbounded')

        if status == gb.GRB.OPTIMAL:
            # get system parameters
            Msys, Dsys, Rsys, Fsys = self._get_sys_papra()
            # get cost
            cost = self.mdl.getObjective().getValue()
            # get fnadir
            nadir = self._get_nadir()
            # get rocof
            rocof = self.dpe / Msys

            print('--------------------- Results -------------------')
            print('Total Cost: %g' % cost)
            print('RoCof prediction: %g ; RoCof limit: %g' % (rocof, self.rocof))
            print('Nadir prediction: %g ; Nadir limit %g' % (nadir, self.fnadir))
    
            # --- get optimization results --
            pg = []
            pru = []
            prd = []
            Mvsg = []
            Dvsg = []
            _, vsg = self._get_GENI_GENII_key()
            for gen in self.gendict.keys():
                pg.append(self.pg[gen].X)
                pru.append(self.pru[gen].X)
                prd.append(self.prd[gen].X)
            for vsg_idx in vsg:
                Mvsg.append(self.Mvsg[vsg_idx].X)
                Dvsg.append(self.Dvsg[vsg_idx].X)
    
        # --- Output dataframe ---
        dcres = pd.DataFrame()
        dcres['gen'] = self.gen['idx']
        dcres['Sn'] = self.gen['Sn']
        dcres['pg'] = pg
        dcres['pru'] = pru
        dcres['prd'] = prd

        MDres = pd.DataFrame()
        MDres['gen'] = vsg
        MDres['Mvsg'] = Mvsg
        MDres['Dvsg'] = Dvsg
        dcres.fillna(0, inplace=True)
        logger.info('Msys and Dsys are normlized by devise Sbase, transform to andes Sbase when do TDS')

        Msys, Dsys, Rsys, Fsys = self._get_sys_papra()
        sys_para = {'Msys': Msys, 'Dsys': Dsys, 'Rsys': Rsys, 'Fsys': Fsys}

        return dcres, MDres, sys_para

    def _get_sys_papra(self):
        """
        Get final system parameters

        Parameters
        ----------
        fnorm : DataFrame
            The DataFrame contains normalization factors for system variables.

        """
        # --- get normalization factors ---
        gendict = self.gendict
        GENI, GENII = self._get_GENI_GENII_key()
        
        # --- Synthetic M/D/F/R ---
        Msys = sum(gendict[gen]['Sn'] * gendict[gen]['M'] for gen in GENI)
        Msys += sum(gendict[gen]['Sn'] * self.Mvsg[gen].X for gen in GENII)
        Msys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Dsys = sum(gendict[gen]['Sn'] * gendict[gen]['D'] for gen in GENI)
        Dsys += sum(gendict[gen]['Sn'] * self.Dvsg[gen].X for gen in GENII)
        Dsys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Rsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * gendict[gen]['Sn'] for gen in GENI)
        Rsys /= sum(gendict[gen]['Sn'] for gen in GENI)

        Fsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * self.pg[gen].X for gen in GENI)
        Fsys /= sum(gendict[gen]['Sn'] for gen in GENI)

        return Msys, Dsys, Rsys, Fsys

    def _get_nadir(self):
        Msys, Dsys, Rsys, Fsys = self._get_sys_papra()
        
        # import network and norm parameters
        fnn = torch.jit.load('net_fnadir.pt')
        fnorm = pd.read_csv('fnorm.csv')

        fmean = fnorm.iloc[0].values
        fstd = fnorm.iloc[1].values

        x = [Msys, Dsys, Rsys, Fsys]

        # norm and predict
        x = (x - fmean[0:4]) / fstd[0:4]
        x = torch.tensor(x, dtype=torch.float32)

        # denorm
        nadir = fnn(x).detach().numpy() * fstd[4] + fmean[4]

        nadir *= self.dpe

        return nadir


# ----------------------------- class: vis2 --------------------------------------
class vis2(vis1):

    def __init__(
                    self, 
                    name='vis2', 
                    norm=None, nn=None, 
                    nn_num=64, 
                    dpe=0.01,
                    rocof_lim = 0.01, 
                    nadir_lim=0.01
                ):

        super().__init__(
                            name = name, 
                            norm = norm, 
                            nn = nn, 
                            nn_num = nn_num, 
                            dpe = dpe, 
                            rocof_lim = rocof_lim, 
                            nadir_lim = nadir_lim
                        )

    def _build_cons(self):
        super()._build_cons()
        self._add_vsg_pcons()

    def _add_vsg_pcons(self):

        gendict = self.gendict
        GENI, GENII = self._get_GENI_GENII_key()
        
        # Synthetic M/D/F/R
        Msys = sum(gendict[gen]['Sn'] * gendict[gen]['M'] for gen in GENI)
        Msys += sum(gendict[gen]['Sn'] * self.Mvsg[gen] for gen in GENII)
        Msys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Dsys = sum(gendict[gen]['Sn'] * gendict[gen]['D'] for gen in GENI)
        Dsys += sum(gendict[gen]['Sn'] * self.Dvsg[gen] for gen in GENII)
        Dsys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Rsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * gendict[gen]['Sn'] for gen in GENI)
        Rsys /= sum(gendict[gen]['Sn'] for gen in GENI)

        Fsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * self.pg[gen] for gen in GENI)
        Fsys /= sum(gendict[gen]['Sn'] for gen in GENI)

        # norm M/D/F/R for nn input
        pnorm = self.norm['pnorm']
        Msys_norm = (Msys - pnorm['M'].iloc[0]) / pnorm['M'].iloc[1] # iloc[0] is mean, iloc[1] is std
        Dsys_norm = (Dsys - pnorm['D'].iloc[0]) / pnorm['D'].iloc[1]
        Fsys_norm = (Fsys - pnorm['Fg'].iloc[0]) / pnorm['Fg'].iloc[1]
        Rsys_norm = (Rsys - pnorm['Rg'].iloc[0]) / pnorm['Rg'].iloc[1]

        hdown = -100
        hup = 100

        # add vsg power constraint for each typeII(vsg) gen
        for gen in GENII:
            # build var
            ap = self.mdl.addVars(self.nn_num, name=f'a_{gen}', vtype=gb.GRB.BINARY)
            zp = self.mdl.addVars(self.nn_num, name=f'z_{gen}', vtype=gb.GRB.CONTINUOUS, obj=0, lb=[0]*self.nn_num)

            # norm Mvsg and Dvsg
            Mvsg_norm = (self.Mvsg[gen] - pnorm['Mvsg'].iloc[0]) / pnorm['Mvsg'].iloc[1]
            Dvsg_norm = (self.Dvsg[gen] - pnorm['Dvsg'].iloc[0]) / pnorm['Dvsg'].iloc[1]

            # add constraint
            # i) nn hidden layer
            zp_bar = []
            for i in range(self.nn_num):
                zp_bar_temp = Msys_norm * self.nn['pw1'][0].iloc[i] + \
                              Dsys_norm * self.nn['pw1'][1].iloc[i] + \
                              Fsys_norm * self.nn['pw1'][2].iloc[i] + \
                              Rsys_norm * self.nn['pw1'][3].iloc[i] + \
                              Mvsg_norm * self.nn['pw1'][4].iloc[i] + \
                              Dvsg_norm * self.nn['pw1'][5].iloc[i] + \
                              self.nn['pb1'][0].iloc[i]
                zp_bar.append(zp_bar_temp)

            self.mdl.addConstrs(zp[i] <= zp_bar[i] - hdown*(1-ap[i]) for i in range(self.nn_num))
            self.mdl.addConstrs(zp[i] >= zp_bar[i] for i in range(self.nn_num))
            self.mdl.addConstrs(zp[i] <= hup * ap[i] for i in range(self.nn_num))

            # ii) nn output layer
            pr_norm = sum(zp[i]*self.nn['pw2'][i].iloc[0] for i in range(self.nn_num)) + self.nn['pb2'].iloc[0]

            pr_pred = pr_norm * pnorm['Ppeak'].iloc[1] + pnorm['Ppeak'].iloc[0]  # denorm

            pr_pred *= self.dpe

            # iii) add constraint
            # Note: 
            # both pru and pud are positive without considering sign
            # but pr_red constains sign, so if dep<0, pr_red == - prd
            if self.dpe > 0:
                self.mdl.addConstr(self.pru[gen] == pr_pred, name=f'{gen}_vsg_pru')
                self.mdl.addConstr(self.pg[gen] + pr_pred <= gendict[gen]['pmax'], name=f'{gen}_vsg_genU')
            else:
                self.mdl.addConstr(self.prd[gen] == -pr_pred, name=f'{gen}_vsg_prd')
                self.mdl.addConstr(self.pg[gen] + pr_pred >= gendict[gen]['pmin'], name=f'{gen}_vsg_genD')

            

    def _get_vsgpr(self, Mvsg, Dvsg):
        Msys, Dsys, Rsys, Fsys = self._get_sys_papra()
        
        # import network and norm parameters
        pnn = torch.jit.load('net_Ppeak.pt')
        pnorm = pd.read_csv('pnorm.csv')

        pmean = pnorm.iloc[0].values
        pstd = pnorm.iloc[1].values

        x = [Msys, Dsys, Rsys, Fsys, Mvsg, Dvsg]

        # norm and predict
        x = (x - pmean[0:6]) / pstd[0:6]
        x = torch.tensor(x, dtype=torch.float32)

        # denorm
        pr = pnn(x).detach().numpy() * pstd[6] + pmean[6]

        pr *= self.dpe
        return pr
        
    def get_res(self):
        """
        Get resutlts, can be used after mdl.optimize().

        Returns
        -------
        DataFrame
            The output DataFrame contains setpoints ``pg``

        """
        self.build()
        self.mdl.optimize()

        status = self.mdl.status
        
        if status == gb.GRB.UNBOUNDED:
            print('The model cannot be solved because it is unbounded')

        if status == gb.GRB.OPTIMAL:
            # get system parameters
            Msys, Dsys, Rsys, Fsys = self._get_sys_papra()
            # get cost
            cost = self.mdl.getObjective().getValue()
            # get fnadir
            nadir = self._get_nadir()
            # get rocof
            rocof = self.dpe / Msys

            print('--------------------- Results -------------------')
            print('Total Cost: %g' % cost)
            print('RoCof prediction: %g ; RoCof limit: %g' % (rocof, self.rocof))
            print('Nadir prediction: %g ; Nadir limit %g' % (nadir, self.fnadir))
    
            # --- get optimization results --
            pg = []
            pru = []
            prd = []
            for gen in self.gendict.keys():
                pg.append(self.pg[gen].X)
                pru.append(self.pru[gen].X)
                prd.append(self.prd[gen].X)

            idx = []
            Sn_vsg = []
            Mvsg = []
            Dvsg = []
            pg_vsg = []
            pru_vsg = []
            prd_vsg = []
            pmax_vsg = []
            pmin_vsg = []
            _, vsg = self._get_GENI_GENII_key()
            for vsg_idx in vsg:
                idx.append(vsg_idx)
                Sn_vsg.append(self.gendict[vsg_idx]['Sn'])
                Mvsg.append(self.Mvsg[vsg_idx].X)
                Dvsg.append(self.Dvsg[vsg_idx].X)
                pg_vsg.append(self.pg[vsg_idx].X)
                pru_vsg.append(self.pru[vsg_idx].X)
                prd_vsg.append(self.prd[vsg_idx].X)
                pmax_vsg.append(self.gendict[vsg_idx]['pmax'])
                pmin_vsg.append(self.gendict[vsg_idx]['pmin'])

        # --- Output results ---
        pgres = pd.DataFrame()
        pgres['gen'] = self.gen['idx']
        pgres['Sn'] = self.gen['Sn']
        pgres['pg'] = pg
        pgres['pru'] = pru
        pgres['prd'] = prd
        pgres.fillna(0, inplace=True)

        vsg_res = pd.DataFrame()
        vsg_res['gen'] = idx
        # vsg_res['Sn'] = Sn_vsg
        vsg_res['Mvsg'] = Mvsg
        vsg_res['Dvsg'] = Dvsg
        vsg_res['pg_vsg'] = pg_vsg
        vsg_res['pru_vsg'] = pru_vsg
        vsg_res['prd_vsg'] = prd_vsg
        vsg_res['pmax_vsg'] = pmax_vsg
        vsg_res['pmin_vsg'] = pmin_vsg
        logger.info('Msys and Dsys are normlized by devise Sbase, transform to andes Sbase when do TDS')

        Msys, Dsys, Rsys, Fsys = self._get_sys_papra()
        sys_para = {'Msys': Msys, 'Dsys': Dsys, 'Rsys': Rsys, 'Fsys': Fsys}

        return pgres, vsg_res, sys_para

# ----------------------- Auxiliary function ----------------------------------

def loadnn(para_path):
    # data path
    dir_path = os.path.abspath('..')
    data_path = dir_path + para_path

    # norm data
    fnorm = pd.read_csv(data_path + '/fnorm.csv')
    pnorm = pd.read_csv(data_path + '/pnorm.csv')
    norm = {'fnorm': fnorm, 'pnorm': pnorm }

    # network parameters
    fw1 = pd.read_csv(data_path + '/fw1.csv', header=None)
    fw2 = pd.read_csv(data_path + '/fw2.csv', header=None)

    fb1 = pd.read_csv(data_path + '/fb1.csv', header=None)
    fb2 = pd.read_csv(data_path + '/fb2.csv', header=None)

    pw1 = pd.read_csv(data_path + '/pw1.csv', header=None)
    pw2 = pd.read_csv(data_path + '/pw2.csv', header=None)

    pb1 = pd.read_csv(data_path + '/pb1.csv', header=None)
    pb2 = pd.read_csv(data_path + '/pb2.csv', header=None)

    nn = {
        'fw1': fw1,
        'fw2': fw2,       
        'fb1': fb1,
        'fb2': fb2,
        'pw1': pw1,
        'pw2': pw2,
        'pb1': pb1,
        'pb2': pb2,
    }
    return nn, norm