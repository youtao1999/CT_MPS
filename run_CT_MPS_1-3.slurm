#!/bin/bash
#SBATCH --job-name=CT_MPS   # Job name with L parameter
#SBATCH --ntasks=4                # Number of tasks (processes)
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --cpus-per-task=1        # Number of CPU cores per task
#SBATCH --time=18:00:00            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ty296/logs/%j.out            # Standard output log
#SBATCH --error=/scratch/ty296/logs/%j.err             # Standard error log
#SBATCH --partition=main


# Variables provided via --export (defaults are handled in Julia script if not provided)
# Disable multithreading to ensure single-threaded workers
export JULIA_NUM_THREADS=1        # Julia threading
export OMP_NUM_THREADS=1          # OpenMP threading (for BLAS/LAPACK)
export MKL_NUM_THREADS=1          # Intel MKL threading
export OPENBLAS_NUM_THREADS=1     # OpenBLAS threading

# Run the command
# singularity exec /scratch/ty296/CT_MPS_mini/julia_CT.sif julia --sysimage /scratch/ty296/CT_MPS_mini/ct_with_wrapper.so /scratch/ty296/CT_MPS_mini/run_CT_MPS_1-3.jl --L $L --p_range "$P_RANGE" --p_fixed_name "$P_FIXED_NAME" --p_fixed_value $P_FIXED_VALUE --ancilla $ANCILLA --maxdim $MAXDIM --threshold $THRESHOLD --n_chunk_realizations $N_CHUNK_REALIZATIONS --job_id $SLURM_JOB_ID --random --output_dir "$OUTPUT_DIR" --store_sv
julia --sysimage=/scratch/ty296/CT_MPS/run_CT_MPS_1-3.so --project=/scratch/ty296/CT_MPS/CT /scratch/ty296/CT_MPS/run_CT_MPS_1-3.jl --L $L --p_ctrl $P_CTRL --p_proj $P_PROJ --seed $SEED 
