using ITensors
using Random
using LinearAlgebra
using MKL
using Pkg
using JSON
using CT
using Printf
using HDF5

using ArgParse
using Serialization
"""
Store a single result directly to HDF5 file, optimized for large singular value arrays.
Uses compression and efficient chunking for large arrays.
Appends results to existing file or creates new file.
"""
function store_result_hdf5_single_shot(filename::String, sv_arr::Vector{Float64}, args::Dict)
    if isfile(filename)
        rm(filename)
    end
    h5open(filename, "cw") do file
        # Handle different data types for chunking
        chunk_size = (min(1000, length(sv_arr)),)
        sv_data = sv_arr
        
        sv_dataset = create_dataset(file, "sv_arr", datatype(Float64), dataspace(sv_data),
        chunk=chunk_size, 
        shuffle=true, deflate=6)
        write(sv_dataset, sv_data)
        
        # Add metadata as attributes to the compressed dataset
        attributes(sv_dataset)["L"] = args["L"]
        attributes(sv_dataset)["p_ctrl"] = args["p_ctrl"]
        attributes(sv_dataset)["p_proj"] = args["p_proj"]
        attributes(sv_dataset)["eps"] = args["eps"]
        attributes(sv_dataset)["seed"] = args["seed"]
    end
end

function main_interactive(L::Int,p_ctrl::Float64,p_proj::Float64,seed::Int, eps::Float64, filename::String)
    args = Dict("L" => L, "p_ctrl" => p_ctrl, "p_proj" => p_proj, "eps" => eps, "seed" => seed)
    ct_f=CT.CT_MPS(L=L,seed=seed,folded=true,store_op=false,store_vec=false,ancilla=0,debug=false,xj=Set([1//3,2//3]),_maxdim=2^(L÷2), _eps=eps, _cutoff=eps)
    i=1
    T_max = 2*(ct_f.L^2)

    for idx in 1:T_max
        println(idx)
        i=CT.random_control!(ct_f,i,p_ctrl,p_proj)
        @show linkdims(ct_f.mps)
        @show "Heap memory usage (MB): " Base.gc_live_bytes() / 1024^2
        @show "Max RSS (MB): " Sys.maxrss() / 1024^2
    end
    sv_arr=CT.von_Neumann_entropy(ct_f.mps,div(ct_f.L,2); n = 0, positivedefinite=false, threshold=1e-16, sv=true)
    store_result_hdf5_single_shot(filename, sv_arr, args)
    return sv_arr, CT.max_bond_dim(ct_f.mps), filename
end

function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--p_ctrl", "-c"
        arg_type = Float64
        default = 0.0
        help = "control rate"
        "--p_proj", "-p"
        arg_type = Float64
        default = 0.0
        help = "measurement rate"
        "--L", "-L"
        arg_type = Int
        default = 8
        help = "system size"
        "--seed", "-s"
        arg_type = Int
        default = 0
        help = "random seed"
        "--eps", "-e"
        arg_type = Float64
        default = 0.0
        help = "set the eps"
        "--output_file", "-o"
        arg_type = String
        default = "test.h5"
        help = "set the output file"
    end
    return parse_args(s)
end

function main()
    println("Uses threads: ",BLAS.get_num_threads())
    println("Uses backends: ",BLAS.get_config())
    args = parse_my_args()
    main_interactive(args["L"], args["p_ctrl"], args["p_proj"], args["seed"],args["eps"], "/scratch/ty296/hdf5_data/p_ctrl0.4_haining/MPS_L$(args["L"])_pctrl$(@sprintf("%.3f", args["p_ctrl"]))_pproj$(@sprintf("%.3f", args["p_proj"]))_s$(args["seed"]).h5")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end