# This script is compute the sum of trajectory fluctuation and the "trajectory-averaged" state fluctuation, given a specific circuit.
import sys
dir_path='../control_transition'
sys.path.append(dir_path)

from tqdm import tqdm
from plot_utils import *
L=10
params_list=[
({'nu':0,'de':1,},
{
# 'p_ctrl':[.47,.49,.51,.53],
# 'p_ctrl':[.6],
'p_ctrl':[0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6],
# 'p_ctrl':[0.1],
'p_proj':np.linspace(0.0,0.0,1),
'sC':np.arange(0,500),
'sm':np.arange(500),
'L':[L]
}
),
]
ob="O"
ob1=ob+'1'
ob2=ob+'2'
for fixed_params,vary_params in params_list:
    data_MPS_0_T_DW_dict=generate_params(
        fixed_params=fixed_params,
        vary_params=vary_params,
        # fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_x01_DW_T.json',
        fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_O_T.json',
        fn_dir_template='/MPS_0-1_C_m_x01_T',
        input_params_template='',
        load_data=load_zip_json,
        filename=None,
        filelist=None,
        load=True,
        data_dict={'fn':set()},
        # zip_fn='/home/jake/Data/MPS_0-1_C_m_x01_T.zip'  
        # zip_fn='./MPS_0-1_C_m_x01_T.zip'  
        # zip_fn=f'./MPS_0-1_C_m_x01_T_L{L}.zip'  
        zip_fn=f'./MPS_0-1_C_m_O_T_L{L}.zip'  
        # zip_fn=f'./MPS_0-1_shots_L{L}.zip' 
    )
df_MPS_0_T_DW=convert_pd(data_MPS_0_T_DW_dict,names=['Metrics','sm','sC','p_ctrl','L','p_proj',])

def traj_qvar(df,L,p_ctrl,sC):
    # trajecotry fluctuation, fix SC, as a function T 
    data = df['observations'].xs(sC,level='sC').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs(L,level='L')
    data_DW1=np.stack(data.xs(ob1,level='Metrics'))
    sigma_mc=data_DW1.var(axis=0)

    # state_fluctuation
    data_DW2=np.stack(data.xs(ob2,level='Metrics'))
    qvar=(data_DW2-data_DW1**2)
    sigma_s=qvar.mean(axis=0)

    return sigma_mc+sigma_s

sigma_shot={}
for p in tqdm(params_list[0][1]['p_ctrl']):
    for L in params_list[0][1]['L']:
        print(p,L)
        sigma_shot[(p,L)]=[]
        for sC in range((params_list[0][1]['sC']).shape[0]):
            # try:
            sigma_shot_=traj_qvar(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC)
            sigma_shot[(p,L)].append(sigma_shot_)
            # except:
            #     pass
        sigma_shot[(p,L)]=np.array(sigma_shot[(p,L)])

# with open(f'traj_state_sum_DW_L{L}.pickle','wb') as f:
with open(f'traj_state_sum_O_L{L}.pickle','wb') as f:
    pickle.dump(sigma_shot,f)