
import os
import andes
import pandas as pd


def get_andes_case(case_path):
    '''
        Get the andes case from the case path.
    
        inputs:
        -------
        case_path: str

        returns:
        --------
        andes attribute
    '''
    # ss0 is used for PP conversion
    dir_path = os.path.abspath('..')

    case = dir_path + case_path
    ssa = andes.load(case,
                    setup=True,
                    no_output=True,
                    default_config=False)

    # Set output mode as 'manual'
    ssa.TDS.config.save_mode = 'manual'

    # set PQ constant load
    ssa.PQ.config.p2p = 1
    ssa.PQ.config.q2q = 1
    ssa.PQ.config.p2z = 0
    ssa.PQ.config.q2z = 0
    ssa.PQ.pq2z = 0

    # Turn on 'numba' to accelerate TDS.
    ssa.config.numba

    return ssa

def get_load(data_path, load_time=10):
    '''
        Get normalized load profile
    
        inputs:
        -------
        data_path: str
        load_time: int
            The time of load profile, default is 10, can also be 18

        returns:
        --------
        load data: pandas dataframe
            normalized by mean value
    '''

    dir_path = os.path.abspath('..')
    path = dir_path + data_path
    d_syn = pd.read_csv(path)

    caseH = load_time
    # the coefficient can be adjusted to fit the case
    if caseH == 10: # load prfile at 10am
        d_syn['sload'] = 1*(d_syn['ha10'] - d_syn['ha10'].min()) / d_syn['ha10'].min() + 1
    if caseH == 18: # load prfile at 6pm
        d_syn['sload'] = 2*(d_syn['ha18'] - d_syn['ha18'].min()) / d_syn['ha18'].min() + 1

    # smooth
    d_syn['sload'] = d_syn['sload'].rolling(20).mean()
        
    # calculate expected load
    step = 300
    d_exp = d_syn.groupby(d_syn.index // step).mean()
    d_exp['time'] = range(0, 3600, 300)

    # align starting point of load with starting point of dispatch results
    d_syn['sload'][0] = d_exp['sload'].iloc[0]
    d_syn['sload'][1:100] = None
    d_syn['sload'] = d_syn['sload'].interpolate(method='polynomial', order=3)

    return d_syn