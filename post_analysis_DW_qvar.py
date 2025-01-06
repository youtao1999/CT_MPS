import sys
dir_path='../control_transition'
sys.path.append(dir_path)
# import matplotlib.pyplot as plt

from tqdm import tqdm
from plot_utils import *

import gc

data_path='.'
# L_list=np.arange(10,41,10)
L=20

# p_ctrl_list=np.round(np.arange(0.4,0.61,0.05),2)
p_ctrl_list= [0.4,0.45,0.47,0.49,0.5,0.51,0.53,0.55,0.6]

params_list=[({'nu':0,'de':1,},{'p_ctrl':p_ctrl_list,'p_proj':np.linspace(0.0,0.0,1),'s':np.arange(10000),'L':[L]}),]

cl_variance_p_ctrl_dict={}

for fixed_params,vary_params in params_list:
    data_MPS_0_DW_dict=generate_params(
        fixed_params=fixed_params,
        vary_params=vary_params,
        fn_template='MPS_({nu},{de})_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_s{s}_x01_DW.json',
        fn_dir_template='./MPS_0-1_DW_x01',
        # fn_dir_template='./MPS_0-1_DW_x12',
        # fn_dir_template='./MPS_0-1_DW_x00',
        input_params_template='',
        load_data=load_zip_json,
        filename=None,
        filelist=None,
        load=True,
        data_dict={'fn':set()},
        zip_fn=os.path.join(data_path,f'MPS_0-1_DW_L{L}.zip')
    )
df_MPS_0_DW=convert_pd(data_MPS_0_DW_dict,names=['Metrics','L','p_ctrl','p_proj','T'])

qvar_dict={}
q_variance_dict={}
q_variance_sem_dict={}

for p_ctrl in tqdm(p_ctrl_list):
    assert np.max(df_MPS_0_DW.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs('DW1',level='Metrics')['observations'].apply(max))<=L, 'issue in the sample'
    qvar=df_MPS_0_DW.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs('DW2',level='Metrics') - df_MPS_0_DW.xs(L,level='L').xs(p_ctrl,level='p_ctrl').xs(0.0,level='p_proj').xs('DW1',level='Metrics').applymap(lambda x : np.array(x)**2)
    qvar_dict[p_ctrl,L]=np.hstack([qvar.loc[t]['observations'] for t in range(L**2,2*L**2)])

    q_variance_dict[p_ctrl,L]=qvar['observations'].apply(np.mean)
    q_variance_sem_dict[p_ctrl,L]=qvar['observations'].apply(np.std)/np.sqrt(qvar['observations'].apply(len))




with open(f'qvar_L{L}.pickle','wb') as f:
    pickle.dump([q_variance_dict,q_variance_sem_dict,qvar_dict],f)


