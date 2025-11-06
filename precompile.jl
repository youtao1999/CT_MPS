# precompile.jl
using CT, ArgParse
include("run_CT_MPS_1-3.jl")
# do one dummy call so all methods get JIT’ed
main_interactive(20, 0.4, 0.7, 42, 0.0, "test.h5")

# To run the precompilation, use the following command:
# export JULIA_DEPOT_PATH=~/julia_depot
# export TMPDIR="/tmp" # (optional)
# using PackageCompiler; using Pkg; Pkg.activate("CT")
# create_sysimage(
#     [:CT, :ITensors, :ArgParse, :JSON],
#     sysimage_path="CT_MPS_1-3_ios.so",
#     precompile_execution_file="precompile.jl"
#   )


# Then to start from the sysimage, use: 
# julia --sysimage ct_with_wrapper.so --project=.
# > using Pkg
# > Pkg.activate("CT")
# > main_interactive(8, 0.1, 0.2, 1, 1)
#
# Or call the script directly:
# julia --sysimage C_m_T.jl.so \
#       --project=CT \
#       -p 2 \
#       run_CT_MPS_C_m_T_multiproc.jl \
#       --params "0.35,0.0,20,99,0,0.35,0.0,20,99,1"