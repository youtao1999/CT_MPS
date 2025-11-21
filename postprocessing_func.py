import numpy as np
from scipy.stats import linregress
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

## helper functions
def regularize(jagged): 
    max_len = max(len(row) for row in jagged)
    arr = np.full((len(jagged), max_len), 1e-6)  # use np.full instead of np.zeros
    for i, row in enumerate(jagged):
        arr[i, :len(row)] = row
    return arr

## observation functions
def corr_len(sv_arr):
    log_sv_arr = np.clip(np.log10(sv_arr), -6, 0)

    # Do the fit and get statistics in one line
    result = linregress(range(len(log_sv_arr)), log_sv_arr)
    if len(log_sv_arr) == 1:
        cxi = 0
        r_squared = 1
    else:
        # Access what you need
        cxi = -1/result.slope
        r_squared = result.rvalue**2  # Note: square the rvalue
    return cxi

def s0(sv_arr):
    return np.log(len(sv_arr))

## ensemble functions
def binder_cumulant(data):
    """Compute Binder cumulant for array of measurements"""
    data = np.array(data)
    m2 = np.mean(data**2)
    m4 = np.mean(data**4)
    U = 1 - m4 / (3 * m2**2)
    return U

def calculate_mean_and_error(arr):
    """Calculate mean and standard error of the mean."""
    array = np.array(arr)
    mean = np.mean(array)
    sem = np.std(array, ddof=1) / np.sqrt(len(array))
    return mean, sem

def calculate_variance_and_error(arr):
    """Calculate variance and standard error of variance."""
    array = np.array(arr)
    mean = np.mean(array)
    var = np.var(array, ddof=1)

    deviations = array - mean
    fourth_moment = np.mean(deviations**4)
    se_var = (1/len(array)) * (fourth_moment - (len(array)-3)/(len(array)-1) * var**2)
    return var, se_var

## export to dataframe
def output_df(data_dict):
    index = pd.MultiIndex.from_tuples([(key[0], key[1]) for key in data_dict.keys() if key!='fn'], names=['p','L'])
    df = pd.DataFrame({'observations': data_dict.values()}, index=index)
    
    # Calculate mean and standard error using helper function
    mean_error = df['observations'].apply(lambda x: calculate_mean_and_error(x))
    df['estimator'] = mean_error.apply(lambda x: x[0])
    df['standard_error'] = mean_error.apply(lambda x: x[1])
    
    # Calculate variance and standard error of variance using helper function
    var_error = df['observations'].apply(lambda x: calculate_variance_and_error(x))
    df['variance'] = var_error.apply(lambda x: x[0])
    df['standard_error_of_variance'] = var_error.apply(lambda x: x[1])
    
    return df
    
def df_obs_from_sv_arr(data_dict, obs_func, ensemble_func):
    data = []
    
    for key in data_dict.keys():
        L, p_proj = key
        obs_lst = []
        
        for sv_arr in data_dict[key]:
            obs_lst.append(obs_func(sv_arr))
        
        estimator_value = ensemble_func(obs_lst)
        data.append((L, p_proj, obs_lst, estimator_value))
    
    # Create MultiIndex from collected data
    index = pd.MultiIndex.from_tuples([(row[0], row[1]) for row in data], 
                                       names=['p', 'L'])
    
    # Create DataFrame with all data
    df = pd.DataFrame({
        'observations': [row[2] for row in data],
        'estimator': [row[3] for row in data]
    }, index=index)
    
    return df

## plotting functions 
def plot_ensemble_metric(df, x_label, y_label, L_list = None):
    sns.set_style("whitegrid")

    # Get sorted L values from the MultiIndex
    L_values = sorted(df.index.get_level_values('L').unique())
    if L_list is not None:
        L_values = [L_val for L_val in L_values if L_val in L_list]
    # Create sequential blue palette (light to dark)
    blues = sns.color_palette("Blues", n_colors=len(L_values)+2)[2:]  # Skip the lightest colors

    fig, ax = plt.subplots(figsize=(10, 6))

    for i, L_val in enumerate(L_values):
        # Select data for this L value from MultiIndex
        group = df.xs(L_val, level='L')
        
        # Get p_proj values from the remaining index level
        subset = group[::2]
        p_proj_values = subset.index.get_level_values('p')
        
        ax.plot(p_proj_values, subset['estimator'], 
                marker='o', label=f'L={L_val}', 
                color=blues[i], linewidth=0.5, markersize=7,
                markeredgewidth=0.5, markeredgecolor='white')

    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    ax.legend(frameon=True)
    plt.tight_layout()
    plt.show()
    
def hist_ensemble_metric(df, x_label, y_label):
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, (L_val, group) in enumerate(df.groupby('L')):
        subset = group[::2]
        ax.hist(subset['estimator'], bins=np.arange(0, 2, 0.1), label=f'L={L_val}', alpha=0.5)
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    ax.legend(frameon=True)