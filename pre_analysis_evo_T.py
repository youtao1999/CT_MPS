import sys
dir_path='../control_transition'
sys.path.append(dir_path)
# import matplotlib.pyplot as plt

from tqdm import tqdm
from plot_utils import *

L=20
params_list=[
({'nu':0,'de':1,},
{
'p_proj':[0,],
'sC':np.arange(0,100),
'sm':[0],
'L':[L],
'maxdim':[2**13,],
# 'cutoff': [1e-10,1e-9,1e-8],
# 'cutoff': [1e-5,1e-3,1e-1],
'cutoff': [1e-15,1e-10, 1e-5,],
'p_ctrl':np.linspace(0,1,21),
}
),
]

for fixed_params,vary_params in params_list:
    data_MPS_0_T_dict=generate_params(
        fixed_params=fixed_params,
        vary_params=vary_params,
        # fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_x01_evo_.json',
        fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_maxdim{maxdim}_cutoff{cutoff:.1e}.json',
        fn_dir_template=f'../CT_MPS/MPS_0-1_evo_L{L}',
        input_params_template='{p_ctrl:.3f},{p_proj:.3f},{L},{sC},{sm},{maxdim},{cutoff:.1e}',
        load_data=load_json,
        filename=f'params_CT_MPS_0_C_m_T_L{L}_series.txt',
        filelist=None,
        load=False,
        data_dict=None,
        # data_dict_file='xj({nu},{de})_C_m.pickle', 
        # data_dict_file='xj({nu},{de})_C_m.json', 
    )

with open(f'params_CT_MPS_0_C_m_T_L{L}_series.txt','r') as f:
    # lines=f.readlines()
    linewidth=100
    count=0
    total_string = []
    string = ''
    for line in f:
        string = string + ',' + line.strip()
        count+=1
        if count>=linewidth:
            count=0
            total_string.append(string[1:])
            string= ''

with open(f'params_CT_MPS_0_C_m_T_L{L}_series.txt','w') as f:
    f.write('\n'.join(total_string))
