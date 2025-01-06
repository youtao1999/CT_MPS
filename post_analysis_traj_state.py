# This script is trying to tell, out of those trajectory fluctuation and state fluctuation being zero, how much of them are due to the state being at fixed point
import sys
dir_path='../control_transition'
sys.path.append(dir_path)
# import matplotlib.pyplot as plt

from tqdm import tqdm
from plot_utils import *
L=10
params_list=[
({'nu':0,'de':1,},
{
# 'p_ctrl':[.47,.49,.51,.53],
'p_ctrl':[.6],
# 'p_ctrl':[0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6],
# 'p_ctrl':[0.1],
'p_proj':np.linspace(0.0,0.0,1),
'sC':np.arange(0,500),
'sm':np.arange(5),
'L':[L]
}
),
]
err=1e-10
for fixed_params,vary_params in params_list:
    data_MPS_0_T_DW_dict=generate_params(
        fixed_params=fixed_params,
        vary_params=vary_params,
        fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_x01_DW_T.json',
        # fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_O_T.json',
        fn_dir_template='/MPS_0-1_C_m_x01_T',
        input_params_template='',
        load_data=load_zip_json,
        filename=None,
        filelist=None,
        load=True,
        data_dict={'fn':set()},
        # zip_fn='/home/jake/Data/MPS_0-1_C_m_x01_T.zip'  
        # zip_fn='./MPS_0-1_C_m_x01_T.zip'  
        zip_fn=f'/home/jake/Data/MPS_0-1_C_m_x01_T_L{L}.zip'  
        # zip_fn=f'./MPS_0-1_C_m_O_T_L{L}.zip'  
        # zip_fn=f'./MPS_0-1_shots_L{L}.zip' 
    )
df_MPS_0_T_DW=convert_pd(data_MPS_0_T_DW_dict,names=['Metrics','sm','sC','p_ctrl','L','p_proj',])

def trajvar(df,L,p_ctrl,sC):
    data=df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs(sC,level='sC')

    data_DW1=np.stack(df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs(sC,level='sC').xs('DW1',level='Metrics')['observations'])
    # fix SC, as a function T 
    mask_traj_0=(data_DW1.std(axis=0)<err)
    mask_fixed_point_traj=(data_DW1.mean(axis=0)<err)
    return mask_traj_0,mask_fixed_point_traj
def qvar(df, L,p_ctrl,sC):
    data_DW1=np.stack(df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs(sC,level='sC').xs('DW1',level='Metrics')['observations'])
    data_DW2=np.stack(df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs(sC,level='sC').xs('DW2',level='Metrics')['observations'])
    qvar=(data_DW2-data_DW1**2)
    mask_qvar_0=(qvar<err)
    mask_fixed_point_qvar=(data_DW1<err)
    return mask_qvar_0,mask_fixed_point_qvar


mask_traj_0_dict={}
mask_fixed_point_traj_dict={}
mask_qvar_0_dict={}
mask_fixed_point_qvar_dict={}
for p in tqdm(params_list[0][1]['p_ctrl']):
    for L in params_list[0][1]['L']:
        print(p,L)
        mask_traj_0_dict[(p,L)]=[]
        mask_fixed_point_traj_dict[(p,L)]=[]
        mask_qvar_0_dict[(p,L)]=[]
        mask_fixed_point_qvar_dict[(p,L)]=[]
        for sC in range((params_list[0][1]['sC']).shape[0]):
            try:
                mask_traj_0,mask_fixed_point_traj=trajvar(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC)
                mask_qvar_0,mask_fixed_point_qvar=qvar(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC)

                mask_traj_0_dict[(p,L)].append(mask_traj_0)
                mask_fixed_point_traj_dict[(p,L)].append(mask_fixed_point_traj)
                mask_qvar_0_dict[(p,L)].append(mask_qvar_0)
                mask_fixed_point_qvar_dict[(p,L)].append(mask_fixed_point_qvar)
            except:
                pass
        mask_traj_0_dict[(p,L)]=np.array(mask_traj_0_dict[(p,L)])
        mask_fixed_point_traj_dict[(p,L)]=np.array(mask_fixed_point_traj_dict[(p,L)])
        mask_qvar_0_dict[(p,L)]=np.array(mask_qvar_0_dict[(p,L)])
        mask_fixed_point_qvar_dict[(p,L)]=np.array(mask_fixed_point_qvar_dict[(p,L)])

with open(f'traj_state_L{L}.pickle','wb') as f:
    pickle.dump([mask_traj_0_dict,mask_fixed_point_traj_dict,mask_qvar_0_dict,mask_fixed_point_qvar_dict],f)