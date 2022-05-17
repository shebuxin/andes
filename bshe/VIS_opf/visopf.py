# vittual inertia scheduling (vis)
# is inherited from dcopf in opf.py

# base class: dcopf
# vis0: dcopf + organize andes data
# vis1: dcopf + Pvsg
# vis2: dcopf + RoCof and fnadir + Pvsg

from statistics import fmean
from andes.interop.pandapower import to_pandapower
from andes.interop.pandapower import make_GSF, build_group_table
import gurobipy as gb
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger(__name__)

from opf.py import dcopf

# ---------------------------------------------------------

class vis0(dcopf):
    """
    vis0: fixed vsg up and down reserve
    """

    def __init__(self, name='rted', norm):
        """
        fnorm: dict
            normalization dict for function fnaidir and ppeak 
            
            fnorm
            {
                Mvsg: [mean, std],
                Dvsg:
            }
        """
        super().__init__(name)
        self.fmean, self.fstd = norm[0]
        self.fstd = norm[1]
        self.fmean = norm[2]
        self.fmean = norm[3]


    def def_type2(self, gen_idx, M_vsg, D_vsg): # used for vsg gen parameter
        """
        Define type 2 generator based on REGCV1/2

        Parameters
        ----------
        gen_idx : str
            Generator index that will be set to type 2.
        M_vsg: float
            Virtual inertia of inverter
        D_vsg: float
            Damping of inverter

        Example:
        ----------
        gen2 = ['PV_4', 'PV_4']
        M_vsg = [1, 1]
        D_vsg = [1, 1]
        ssd.def_type2(gen2, M_vsg, D_vsg)
        """

        self.gen['type'] = 1
        self.gen['Mvsg'] = 0
        self.gen['Dvsg'] = 0
        for idx, M, D in zip(gen_idx, M_vsg, D_vsg):
            row = self.gen[self.gen['idx'] == idx].index[0]
            self.gen['type'].iloc[row] = 2
            self.gen['Mvsg'].iloc[row] = M
            self.gen['dvsg'].iloc[row] = D
    

    def build(self):
        self.data_check()

        # --- build RTED model ---
        self.update_dict()
        self.mdl = gb.Model(self.name)
        self.mdl = self._build_vars(self.mdl)
        self.mdl = self._build_obj(self.mdl)
        self.mdl = self._build_cons(self.mdl)
        logger.info('Successfully build vis0 model.')
