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
# 'p_ctrl':[.5,],
'p_ctrl':[0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6],
'p_proj':np.linspace(0.0,0.0,1),
# 'sC':np.arange(0,500),
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
        # fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_O_T.json',
        fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_sC{sC}_sm{sm}_shots.json',
        fn_dir_template='/MPS_0-1_C_m_T',
        input_params_template='',
        load_data=load_zip_json,
        filename=None,
        filelist=None,
        load=True,
        data_dict={'fn':set()},
        # zip_fn=f'./MPS_0-1_C_m_T_L{L}.zip'  
        # zip_fn=f'./MPS_0-1_C_m_O_T_L{L}.zip' 
        zip_fn=f'./MPS_0-1_shots_L{L}.zip' 
    )
df_MPS_0_T_DW=convert_pd(data_MPS_0_T_DW_dict,names=['Metrics','sm','sC','p_ctrl','L','p_proj',])
def circ_var_dw(df,L,p_ctrl,sC):
    # try:
    data=df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0,level='p_proj').xs(sC,level='sC')
    # return np.stack(data.xs('DW1')['observations'].values).mean(axis=0)
    return np.stack(data.xs('O1')['observations'].values).mean(axis=0)
    # except:
        # print(f'Missing data for L={L}, p_ctrl={p_ctrl}, sC={sC}')
        # return None
def circ_var_dw_all(df,L,p_ctrl,sC):
    try:
        return np.vstack(df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0,level='p_proj').xs('O1',level='Metrics').xs(sC,level='sC')['observations']).mean(axis=0)
    except:
        print(f'Missing data for L={L}, p_ctrl={p_ctrl}, sC={sC}')
        return None


circ_var_dw_dict={}
circ_var_dw_sem_dict={}
for p in tqdm(params_list[0][1]['p_ctrl']):
    print(p,L)
    sC_circ_var_dw=[circ_var_dw_all(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC) for sC in np.arange((params_list[0][1]['sC']).shape[0])]
    sC_circ_var_dw=np.array([x for x in sC_circ_var_dw if x is not None])
    sC_circ_var=np.var(sC_circ_var_dw,axis=0)
    sC_circ_var_sem=sC_circ_var*np.sqrt(2/((params_list[0][1]['sC']).shape[0]-1))
    circ_var_dw_dict[p,L],circ_var_dw_sem_dict[p,L] = sC_circ_var, sC_circ_var_sem

# with open(f'circ_var_C_m_T_L{L}.pickle','wb') as f:
# with open(f'circ_var_C_m_T_O_L{L}.pickle','wb') as f:
with open(f'circ_var_shots_L{L}.pickle','wb') as f:
    pickle.dump([circ_var_dw_dict,circ_var_dw_sem_dict,],f)