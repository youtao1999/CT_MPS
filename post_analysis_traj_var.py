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
# 'p_ctrl':[.47,.49,.51,.53],
# 'p_ctrl':[.4,.5,.6],
# 'p_ctrl':[0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6],
# 'p_ctrl':[0.1],
'p_ctrl':[0.],
# 'p_proj':np.linspace(0.0,0.0,1),
'p_proj':np.round(np.arange(0.4,1,0.1),2),
'sC':np.arange(0,500),
'sm':np.arange(500),
'L':[L]
}
),
]
ob="DW"
ob1=ob+'1'
ob2=ob+'2'
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
        # zip_fn=f'./MPS_0-1_C_m_x01_T_L{L}.zip'
        # zip_fn=f'/home/jake/Data/MPS_0-1_C_m_x01_T_L{L}.zip'  
        # zip_fn=f'./MPS_0-1_C_m_O_T_L{L}.zip'  
        # zip_fn=f'./MPS_0-1_shots_L{L}.zip' 
        zip_fn=f'./MPS_0-1_MIPT_T_L{L}.zip' 
    )
df_MPS_0_T_DW=convert_pd(data_MPS_0_T_DW_dict,names=['Metrics','sm','sC','p_ctrl','L','p_proj',])

# def trajvar(df,L,p_ctrl,sC):
#     data=df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs(sC,level='sC')

#     # single=[data.xs(sm,level='sm').loc['DW1']['observations'] for sm in range((params_list[0][1]['sm']).shape[0])]
#     single=[data.xs(sm,level='sm').loc[ob1]['observations'] for sm in range((params_list[0][1]['sm']).shape[0])]

#     return np.array(single).var(axis=0)

def trajvar(df,L,p_ctrl,sC,p_proj=0):
    data = df['observations'].xs(sC,level='sC').xs(p_ctrl,level='p_ctrl').xs(p_proj,level='p_proj').xs(L,level='L')
    data_DW1=np.stack(data.xs(ob1,level='Metrics'))
    sigma_mc=data_DW1.var(axis=0)
    return sigma_mc

traj_var_dw_dict={}
traj_var_dw_sem_dict={}
sC_traj_var_dw={}
# for p in tqdm(params_list[0][1]['p_ctrl']):
for p in tqdm(params_list[0][1]['p_proj']):
    for L in params_list[0][1]['L']:
        print(p,L)
        # sC_traj_var_dw=[trajvar(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC) for sC in range((params_list[0][1]['sC']).shape[0])]
        sC_traj_var_dw[(p,L)]=[]
        for sC in range((params_list[0][1]['sC']).shape[0]):
            try:
                sC_traj_var_dw[(p,L)].append(trajvar(df_MPS_0_T_DW,L=L,p_ctrl=0,sC=sC,p_proj=p))
            except:
                pass
        sC_traj_var_dw[(p,L)]=np.array(sC_traj_var_dw[(p,L)])
        traj_var_dw_dict[(p,L)], traj_var_dw_sem_dict[(p,L)]=np.mean(sC_traj_var_dw[(p,L)],axis=0), np.std(sC_traj_var_dw[(p,L)],axis=0)/np.sqrt((params_list[0][1]['sC']).shape[0])

# with open(f'C_m_T_L{L}.pickle','wb') as f:
# with open(f'traj_var_C_m_T_O_L{L}.pickle','wb') as f:
# with open(f'circ_var_shots_L{L}.pickle','wb') as f:
with open(f'traj_var_C_m_T_MIPT_L{L}.pickle','wb') as f:
    pickle.dump([sC_traj_var_dw,traj_var_dw_dict,traj_var_dw_sem_dict],f)