import sys
dir_path='../control_transition'
sys.path.append(dir_path)
# import matplotlib.pyplot as plt

from tqdm import tqdm
from plot_utils import *

L=40
params_list=[
({'nu':0,'de':1,},
{
# 'p_ctrl':[.47,.49,.51,.53],
# 'p_ctrl':[.4,.5,.6],
'p_ctrl':[0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6],
'p_proj':np.linspace(0.0,0.0,1),
'sC':np.arange(0,500),
'sm':np.arange(500),
'L':[L,]
# 'L':[40,]
}
),
]

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
        # zip_fn='./MPS_0-1_C_m_x01_T_L40.zip'  
        zip_fn=f'./MPS_0-1_C_m_O_T_L{L}.zip'  
    )
df_MPS_0_T_DW=convert_pd(data_MPS_0_T_DW_dict,names=['Metrics','sm','sC','p_ctrl','L','p_proj',])

def qvar_dw(df, L,p_ctrl,sC):
    data=df.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0,level='p_proj').xs(sC,level='sC')
    # single=np.array([data.xs(sm,level='sm').loc['DW2']['observations']-data.xs(sm,level='sm').loc['DW1']['observations']**2 for sm in range((params_list[0][1]['sm']).shape[0])])
    single=np.array([data.xs(sm,level='sm').loc['O2']['observations']-data.xs(sm,level='sm').loc['O1']['observations']**2 for sm in range((params_list[0][1]['sm']).shape[0])])
    return single.mean(axis=0)

qvar_dw_dict={}
qvar_dw_sem_dict={}
sC_qvar_dw_dict={}
for p in tqdm(params_list[0][1]['p_ctrl']):
    for L in params_list[0][1]['L']:
        print(p,L)
        # sC_qvar_dw=[qvar_dw(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC) for sC in np.arange((params_list[0][1]['sC']).shape[0])]
        sC_qvar_dw_dict[(p,L)]=[]
        for sC in range((params_list[0][1]['sC']).shape[0]):
            try:
                sC_qvar_dw_dict[(p,L)].append(qvar_dw(df_MPS_0_T_DW,L=L,p_ctrl=p,sC=sC))
            except:
                pass
        sC_qvar_dw_dict[(p,L)]=np.array(sC_qvar_dw_dict[(p,L)])
        qvar_dw_dict[(p,L)]=np.mean(sC_qvar_dw_dict[(p,L)],axis=0)
        qvar_dw_sem_dict[(p,L)]=np.std(sC_qvar_dw_dict[(p,L)],axis=0)/np.sqrt((params_list[0][1]['sC']).shape[0]*(params_list[0][1]['sm']).shape[0])
        
# with open('qvar_C_m_T_L30.pickle','wb') as f:
with open(f'qvar_C_m_T_O_L{L}.pickle','wb') as f:
    pickle.dump([qvar_dw_dict,qvar_dw_sem_dict,sC_qvar_dw_dict],f)