using Distributed

# Add worker processes (do this early, before loading other packages)
# Using 5 total processes: 1 main + 4 workers
addprocs(4)

# Now load packages on all workers
@everywhere using ITensors
@everywhere using Random
@everywhere using LinearAlgebra
@everywhere using MKL
@everywhere using JSON
@everywhere using CT
@everywhere using Printf
@everywhere using HDF5
@everywhere using ArgParse
@everywhere using Serialization

# Define your functions with @everywhere so they're available on all workers
@everywhere function store_result_hdf5_single_shot(filename::String, sv_arr::Vector{Float64}, args::Dict)
    if isfile(filename)
        rm(filename)
    end
    h5open(filename, "cw") do file
        chunk_size = (min(1000, length(sv_arr)),)
        sv_data = sv_arr
        sv_dataset = create_dataset(file, "sv_arr", datatype(Float64), dataspace(sv_data),
            chunk=chunk_size, 
            shuffle=true, deflate=6)
        write(sv_dataset, sv_data)
        
        attributes(sv_dataset)["L"] = args["L"]
        attributes(sv_dataset)["p_ctrl"] = args["p_ctrl"]
        attributes(sv_dataset)["p_proj"] = args["p_proj"]
        attributes(sv_dataset)["eps"] = args["eps"]
        attributes(sv_dataset)["seed"] = args["seed"]
    end
end

@everywhere function main_interactive(L::Int, p_ctrl::Float64, p_proj::Float64, seed::Int, eps::Float64, filename::String)
    args = Dict("L" => L, "p_ctrl" => p_ctrl, "p_proj" => p_proj, "eps" => eps, "seed" => seed)
    ct_f = CT.CT_MPS(L=L, seed=seed, folded=true, store_op=false, store_vec=false, 
                     ancilla=0, debug=false, xj=Set([1//3, 2//3]), 
                     _maxdim=2^(L÷2), _eps=eps, _cutoff=eps)
    i = 1
    T_max = 2*(ct_f.L^2)
    
    for idx in 1:T_max
        println("Worker $(myid()): seed $seed, step $idx")
        i = CT.random_control!(ct_f, i, p_ctrl, p_proj)
    end
    
    sv_arr = CT.von_Neumann_entropy(ct_f.mps, div(ct_f.L, 2); 
                                     n=0, positivedefinite=false, 
                                     threshold=1e-16, sv=true, cutoff=eps)
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
        "--seeds", "-s"
            arg_type = Int
            nargs = '*'
            default = [0, 1, 2, 3, 4]
            help = "list of random seeds to run"
        "--eps", "-e"
            arg_type = Float64
            default = 1e-10
            help = "set the eps"
        "--output_dir", "-o"
            arg_type = String
            default = "/scratch/ty296/hdf5_data/p_ctrl0.4_haining/cutoff1e-10/"
            help = "output directory"
    end
    return parse_args(s)
end

function main()
    println("Main process uses threads: ", BLAS.get_num_threads())
    println("Main process uses backends: ", BLAS.get_config())
    println("Total processes: ", nprocs(), " (1 main + ", nworkers(), " workers)")
    
    args = parse_my_args()
    
    # Run multiple simulations in parallel with different seeds
    seeds = args["seeds"]
    
    results = @distributed (vcat) for seed in seeds
        filename = joinpath(args["output_dir"], 
                           "MPS_L$(args["L"])_pctrl$(@sprintf("%.3f", args["p_ctrl"]))_pproj$(@sprintf("%.3f", args["p_proj"]))_s$(seed).h5")
        
        sv_arr, max_bond, file = main_interactive(args["L"], args["p_ctrl"], 
                                                   args["p_proj"], seed, 
                                                   args["eps"], filename)
        
        # Return summary info
        [(seed=seed, max_bond=max_bond, filename=file)]
    end
    
    println("\nAll simulations completed!")
    println("Results summary:")
    for r in results
        println("  Seed $(r.seed): max_bond=$(r.max_bond), file=$(r.filename)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# julia --project=CT --sysimage=run_CT_MPS_1-3.so run_CT_MPS_1-3_dist.jl --L 12 --p_ctrl 0.4 --p_proj 0.7 --seeds 123 124 125 126 --output_dir /scratch/ty296/test_output/