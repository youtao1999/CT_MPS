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
# 'p_ctrl':[.5,],
# 'p_ctrl':[0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6],
'p_ctrl':[.65,.7,.75,.8,.85],
'p_proj':np.linspace(0.0,0.0,1),
# 'sC':np.arange(0,5),
'sC':np.arange(0,500),
'sm':np.arange(500),
'L':[L]
# 'L':[40,]
}
),
]

for fixed_params,vary_params in params_list:
    data_MPS_0_T_DW_dict=generate_params(
        fixed_params=fixed_params,
        vary_params=vary_params,
        # fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_DW_T.json',
        fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_O_T.json',
        fn_dir_template='/MPS_0-1_C_m_T',
        input_params_template='',
        load_data=load_zip_json,
        filename=None,
        filelist=None,
        load=True,
        data_dict={'fn':set()},
        # zip_fn=f'./MPS_0-1_C_m_T_L{L}.zip'  
        zip_fn=f'./MPS_0-1_C_m_O_T_L{L}.zip' 
    )
df_MPS_0_T_DW=convert_pd(data_MPS_0_T_DW_dict,names=['Metrics','sm','sC','p_ctrl','L','p_proj',])
def circ_var_dw_all(df,L,p_ctrl,sC):
    try:
        return np.vstack(df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0,level='p_proj').xs('O1',level='Metrics').xs(sC,level='sC')['observations'])
    except:
        print(f'Missing data for L={L}, p_ctrl={p_ctrl}, sC={sC}')
        return None


O_mean_dict={}
O_sem_dict={}
for p in tqdm(params_list[0][1]['p_ctrl']):
    print(p,L)
    sC_all=[circ_var_dw_all(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC) for sC in np.arange((params_list[0][1]['sC']).shape[0])]
    data = np.vstack([x for x in sC_all if x is not None])
    

    O_mean_dict[p,L]= data.mean(axis=0)
    O_sem_dict[p,L]= data.std(axis=0)/np.sqrt(data.shape[0])

with open(f'mean_C_m_T_O_L{L}.pickle','wb') as f:
    pickle.dump([O_mean_dict,O_sem_dict],f)