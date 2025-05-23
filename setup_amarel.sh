#!/bin/bash

# Define variables
AMAREL_USER="ty296"  # Your Amarel username
AMAREL_HOST="amarel.rutgers.edu"
REMOTE_SCRATCH="/scratch/${AMAREL_USER}"
CONTAINER_NAME="julia_itensors.sif"

# Check if we're doing local testing
if [ "$1" == "--local" ]; then
    echo "Building container locally..."
    singularity build --fakeroot ${CONTAINER_NAME} julia_itensors.def
    echo "Local build complete! Container created as ${CONTAINER_NAME}"
    exit 0
fi

# Step 1: Create remote directories if they don't exist
ssh amarel "mkdir -p ${REMOTE_SCRATCH}/CT_MPS"

# Step 2: Copy source code and definition files to Amarel
echo "Copying files to Amarel..."
rsync -av --progress \
    ./CT \
    ./run_CT_MPS.jl \
    ./run_CT_MPS_coherence.jl \
    ./julia_itensors.def \
    ${AMAREL_USER}@${AMAREL_HOST}:${REMOTE_SCRATCH}/CT_MPS/

# Step 3: Build Singularity container on Amarel using fakeroot
echo "Building container on Amarel (this may take a while)..."
ssh amarel "cd ${REMOTE_SCRATCH}/CT_MPS && \
    module load singularity && \
    singularity build --fakeroot ${CONTAINER_NAME} julia_itensors.def"

# Step 4: Create a test script
ssh amarel "cd ${REMOTE_SCRATCH}/CT_MPS && cat > test_run.sh << 'EOL'
#!/bin/bash
#SBATCH --job-name=CT_MPS_test
#SBATCH --output=CT_MPS_test_%j.out
#SBATCH --error=CT_MPS_test_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=16G

module load singularity

cd ${REMOTE_SCRATCH}/CT_MPS
singularity exec ${CONTAINER_NAME} julia run_CT_MPS.jl --p_ctrl 0.5 --p_proj 0.0 --L 8 --seed 0 --ancilla 0
EOL"

# Make the test script executable
ssh amarel "chmod +x ${REMOTE_SCRATCH}/CT_MPS/test_run.sh"

echo "Setup complete! To run a test job, connect to Amarel and execute:"
echo "cd ${REMOTE_SCRATCH}/CT_MPS && sbatch test_run.sh" 