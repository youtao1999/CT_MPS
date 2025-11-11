#!/bin/bash
#SBATCH --job-name=CT_MPS   # Job name with L parameter
#SBATCH --ntasks=1                # Number of tasks (processes)
#SBATCH --requeue
#SBATCH --cpus-per-task=5        # Number of CPU cores per task (1 main + 4 workers)
#SBATCH --mem-per-cpu=16G
#SBATCH --time=36:00:00            # Time limit hrs:min:sec
#SBATCH --output=/scratch/ty296/logs/%j.out            # Standard output log
#SBATCH --error=/scratch/ty296/logs/%j.err             # Standard error log
#SBATCH --partition=main


# Disable multithreading to ensure single-threaded workers
export JULIA_NUM_THREADS=1        # Julia threading
export OMP_NUM_THREADS=1          # OpenMP threading (for BLAS/LAPACK)
export MKL_NUM_THREADS=1          # Intel MKL threading
export OPENBLAS_NUM_THREADS=1     # OpenBLAS threading

# Variables provided via --export (defaults are handled in Julia script if not provided)

# Run the command
# singularity exec /scratch/ty296/CT_MPS_mini/julia_CT.sif julia --sysimage /scratch/ty296/CT_MPS_mini/ct_with_wrapper.so /scratch/ty296/CT_MPS_mini/run_CT_MPS_1-3.jl --L $L --p_range "$P_RANGE" --p_fixed_name "$P_FIXED_NAME" --p_fixed_value $P_FIXED_VALUE --ancilla $ANCILLA --maxdim $MAXDIM --threshold $THRESHOLD --n_chunk_realizations $N_CHUNK_REALIZATIONS --job_id $SLURM_JOB_ID --random --output_dir "$OUTPUT_DIR" --store_sv
julia --sysimage=/scratch/ty296/CT_MPS/run_CT_MPS_1-3.so --project=/scratch/ty296/CT_MPS/CT /scratch/ty296/CT_MPS/run_CT_MPS_1-3_dist.jl --L $L --p_ctrl $P_CTRL --p_proj $P_PROJ --seeds $SEEDS --output_dir $OUTPUT_DIR --eps $EPS

# sbatch --export=L=12,P_CTRL=0.4,P_PROJ=0.7,SEEDS="1 2 3 4 5 6 7 8 9 10",OUTPUT_DIR=/scratch/ty296/test_output/ /scratch/ty296/CT_MPS/run_CT_MPS_1-3_dist.slurm
