#!/bin/bash
#SBATCH --job-name=CT_MPS_ensemble
#SBATCH --output=CT_MPS_ensemble_%A.out
#SBATCH --error=CT_MPS_ensemble_%A.err
#SBATCH --ntasks=200  # Default number of tasks
#SBATCH --cpus-per-task=1  # Ensure one CPU per task
#SBATCH --time=24:00:00
#SBATCH --mem=64G

# Default values if not provided via --export
: ${L:=8}
: ${SCAN_TYPE:="p_ctrl"}
: ${P_FIXED:=0.0}
: ${P_RANGE:="0.45:0.55:10"}
: ${ANCILLA:=0}
: ${N_REAL:=2000}

# Load any required modules (modify as needed for your cluster)
module purge
module load singularity

# Set environment variables to ensure single-threading
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export JULIA_WORKER_TIMEOUT=3600  # Increase worker timeout to 1 hour

# Run the Julia script with distributed processing
# The -p flag tells Julia how many worker processes to start (should match --ntasks)
singularity exec julia_itensor.sif julia --check-bounds=no -O3 -p $SLURM_NTASKS run_CT_MPS_ensemble.jl \
    --scan_type "$SCAN_TYPE" \
    --p_fixed $P_FIXED \
    --p_range "$P_RANGE" \
    --L $L \
    --ancilla $ANCILLA \
    --n_realizations $N_REAL 


# Example usage: sbatch --export=L=8,NTASKS=100,P_RANGE="0.1:0.9:20" run_ensemble.slurm

# for interactive session
# first request resource allocation
# srun --nodes=1 --ntasks=10 --cpus-per-task=1 --mem=40G --time=01:00:00 --pty bash
# make sure to load singularity module
# module load singularity
# then run the test
# singularity exec julia_itensor.sif julia --threads=1 -p 10 run_CT_MPS_ensemble.jl --scan_type p_ctrl --p_fixed 0.0 --p_range "0.45:0.55:3" --L 8 --ancilla 0 --n_realizations 100 --maxdim 10