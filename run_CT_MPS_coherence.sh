#!/bin/bash
echo "Hello OSPool from Job running on `whoami`@`hostname`" 
mkdir -p CT/src
mv Project.toml CT
mv Manifest.toml CT
mv CT.jl CT/src

# export JULIA_DEPOT_PATH=/usr/local/julia/share/julia
export JULIA_DEPOT_PATH=/srv/.julia:/usr/local/julia/share/julia

# julia --sysimage /usr/local/sysimage.so run_CT_MPS.jl --p $1 --L $2 --seed $3 --ancilla $4
julia run_CT_MPS_coherence.jl --p_ctrl $1 --p_proj $2 --L $3 --seed $4
# julia -e "using JSON; println(\"Success loading JSON\")"

# julia  run_CT_MPS.jl --p 1. --L 8 --seed 1 --ancilla 1
# julia --sysimage /usr/local/julia/share/julia/sysimages/sys_itensors.so run_CT_MPS.jl --p 1. --L 8 --seed 1 --ancilla 1