#!/bin/bash
echo "Hello OSPool from Job running on `whoami`@`hostname`" 
mkdir -p CT/src
mv Project.toml CT
mv Manifest.toml CT
mv CT.jl CT/src

export JULIA_DEPOT_PATH=/srv/.julia:/usr/local/julia/share/julia

julia run_CT_MPS_C_m.jl --p_ctrl $1 --p_proj $2 --L $3 --seed_C $4 --seed_m $5