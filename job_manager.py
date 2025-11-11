# python script to monitor jobs statues in queue

from typing import Any


import subprocess
import time
import tempfile
import json
import os
import numpy as np
from itertools import product
from collections import defaultdict


if __name__ == "__main__":

    # Configuration parameters
    
    L_list = [12, 16, 20, 24]
    p_ctrl_list = [0.4]
    p_proj_list = np.linspace(0.5, 1.0, 50)
    seed_list = range(2000, 4000, 1)
    all_combinations = list(product(L_list, p_ctrl_list, p_proj_list, seed_list))
    
    # filter out existing files
    job_list = []
    for L, p_ctrl, p_proj, seed in all_combinations:
        file_name = f"/scratch/ty296/hdf5_data/p_ctrl0.4_haining/MPS_L{L}_pctrl{p_ctrl:.3f}_pproj{p_proj:.3f}_s{seed}.h5"
        if os.path.exists(file_name):
            print(f"File {file_name} exists, skipping", flush=True)
        else:
            job_list.append((L, p_ctrl, p_proj, seed))
    
    total_simulations = len(job_list)
    print(f"Jobs to submit: {total_simulations} jobs", flush=True)

    max_queue_size = 480
    time_interval = 60
    SLURM_SCRIPT = "/scratch/ty296/CT_MPS/run_CT_MPS_1-3.slurm"
    OUTPUT_DIR = "/scratch/ty296/hdf5_data/p_ctrl0.4_haining/cutoff1e-10/"

    total_num_jobs = len(job_list)
    # submit jobs
    while True:
        job_ids = subprocess.check_output("squeue -u ty296 --format=%j", shell=True).decode("utf-8").split("\n")
        job_ids = [job_id.strip() for job_id in job_ids if job_id.strip()]
        print(f'"jobs in queue: {len(job_ids)}"', flush=True)
        if len(job_ids) < max_queue_size:
            if len(job_list) == 0:
                break
            else:
                num_jobs_to_submit = min(max_queue_size - len(job_ids), len(job_list))
                for job in job_list[:num_jobs_to_submit]:
                    L, p_ctrl, p_proj, seed = job
                    cmd = f'sbatch --export=ALL,L={L},P_CTRL={p_ctrl},P_PROJ={p_proj},SEED={seed},OUTPUT_DIR={OUTPUT_DIR} {SLURM_SCRIPT}'
                    print(f"Command: {cmd}", flush=True)
                    subprocess.run(cmd, shell=True)
                print(f'Submitted {num_jobs_to_submit} jobs to queue', flush=True)
                # Remove submitted jobs from the list
                job_list = job_list[num_jobs_to_submit:]
                # Show progress
                percentage_complete = (total_simulations - len(job_list)) / total_simulations * 100
                print(f'Progress: {total_simulations - len(job_list)}/{total_simulations} simulations ({percentage_complete:.2f}% complete, {len(job_list)} jobs remaining)', flush=True)
        else:
            time.sleep(time_interval)

# nohup python3 -u /scratch/ty296/CT_MPS/job_manager.py > /scratch/ty296/CT_MPS/job_manager.log 2>&1 &

