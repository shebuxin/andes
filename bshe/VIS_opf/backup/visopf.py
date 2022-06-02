from os import PRIO_PGRP
from statistics import fmean
from andes.interop.pandapower import to_pandapower
from andes.interop.pandapower import make_GSF, build_group_table
import gurobipy as gb
import pandas as pd
import numpy as np
import logging
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
'''

# ---------------------------------------------------------

class vis1(dcopf):
    """
    vis1: dcopf + fnadir/RoCof (ML linearization)
    """

    def __init__(self, name='vis1', norm=None, nn=None, nn_num=64, dpe=0.01):
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
        super().__init__(name)
        self.norm = norm 
        self.nn = nn
        self.nn_num = nn_num
        self.fnadir = 0.01 # 0.6Hz
        self.rocof = 0.01 # 0.6Hz/s
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

        self.Msys_test = Msys

        Dsys = sum(gendict[gen]['Sn'] * gendict[gen]['D'] for gen in GENI)
        Dsys += sum(gendict[gen]['Sn'] * self.Dvsg[gen] for gen in GENII)
        Dsys /= sum(gendict[gen]['Sn'] for gen in gendict.keys())

        Rsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * gendict[gen]['Sn'] for gen in GENI)
        Rsys /= sum(gendict[gen]['Sn'] for gen in GENI)

        Fsys = sum(gendict[gen]['K'] / gendict[gen]['R'] * self.pg[gen] for gen in GENI)
        Fsys /= sum(gendict[gen]['Sn'] for gen in GENI)
        
        # --- add gurobi var for fnadir ---
        af = []
        zf = []
        for i in range(self.nn_num):
            af.append('af'+str(i))
            zf.append('zf'+str(i))
        self.af_key = af
        self.zf_key = zf
        self.af = self.mdl.addVars(af, name='af', vtype=gb.GRB.BINARY)
        self.zf = self.mdl.addVars(zf, name='zf', vtype=gb.GRB.CONTINUOUS, obj=0,
                                lb=[0]*self.nn_num)

        # --- add constraints ---
        # --- RoCof con----
        self.mdl.addConstr(self.rocof * Msys >= abs(self.dpe), name='RoCof')

        # --- fnadir con----
        # - self.fnadir <= fnadir <= self.fnadir
        fnorm = self.norm['fnorm']
        Msys_norm = (Msys - fnorm['M'].iloc[0]) / fnorm['M'].iloc[1]
        Dsys_norm = (Dsys - fnorm['D'].iloc[0]) / fnorm['D'].iloc[1]
        Fsys_norm = (Fsys - fnorm['Fg'].iloc[0]) / fnorm['Fg'].iloc[1]
        Rsys_norm = (Rsys - fnorm['Rg'].iloc[0]) / fnorm['Rg'].iloc[1]

        zf_bar = []
        num_idx = []
        for i in range(self.nn_num):
            zf_bar_temp = Msys_norm * self.nn['fw1'][0].iloc[i] + \
                          Dsys_norm * self.nn['fw1'][1].iloc[i] + \
                          Fsys_norm * self.nn['fw1'][2].iloc[i] + \
                          Rsys_norm * self.nn['fw1'][3].iloc[i] + self.nn['fb1'][0].iloc[i]
            zf_bar.append(zf_bar_temp)
            num_idx.append(i)
        
        self.zf_bar = zf_bar

        hdown = -100
        hup = 100

        for zf_key, af_key, zf_bar_ in zip(self.zf_key, self.af_key, zf_bar):
            self.mdl.addConstr( self.zf[zf_key] <= zf_bar_ - hdown * (1 - self.af[af_key]) )
            self.mdl.addConstr( self.zf[zf_key] >= zf_bar_ )
            self.mdl.addConstr( self.zf[zf_key] <= hup * self.af[af_key] )

        fnadir_norm = 0.0
        for zf_key, i in zip(self.zf_key, num_idx):
            fnadir_norm += self.zf[zf_key] * self.nn['fw2'][i].iloc[0]
        fnadir_norm += self.nn['fb2'][0].iloc[0]

        fnadir_pred = fnadir_norm * fnorm['fnadir'].iloc[1] + fnorm['fnadir'].iloc[0]
        fnadir_pred *= self.dpe

        self.mdl.addConstr(fnadir_pred >= - self.fnadir, name=f'fnadir_D')
        self.mdl.addConstr(fnadir_pred <= self.fnadir, name=f'fnadir_U')

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
        # --- check if mdl is sovled ---
        if not hasattr(self.pg[self.gen.idx[0]], 'X'):
            logger.warning('vis1 has no valid resutls!')
            pg = [0] * self.gen.shape[0]
            pru = [0] * self.gen.shape[0]
            prd = [0] * self.gen.shape[0]
        else:
            logger.warning('Successfully solve vis1.')

            # --- gather data --
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
            # --- cost ---
            self.res_cost = self.mdl.getObjective().getValue()
            logger.info(f'Total cost={np.round(self.res_cost, 3)}')
            
        # --- build output table ---
        dcres = pd.DataFrame()
        dcres['gen'] = self.gen['idx']
        dcres['pg'] = pg
        dcres['pru'] = pru
        dcres['prd'] = prd

        MDres = pd.DataFrame()
        MDres['gen'] = vsg
        MDres['Mvsg'] = Mvsg
        MDres['Dvsg'] = Dvsg
        dcres.fillna(0, inplace=True)
        logger.info('Msys and Dsys are normlized by devise Sbase, transform to andes Sbase when do TDS')

        return dcres, MDres