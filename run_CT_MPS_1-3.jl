using Pkg
Pkg.activate("CT")
using CT
using ITensors
using ITensorMPS
using Random
using LinearAlgebra
using MKL
using Pkg
using JSON
using Printf
using ProgressMeter

using ArgParse
using Serialization

function get_mem_stats()
    stats = Base.gc_num()
    return Dict(
        "allocated_bytes" => stats.total_allocd,
        "gc_total_time" => stats.total_time/1e9,  # Convert to seconds
        "gc_pause_time" => stats.pause/1e9,  # Convert to seconds
        "gc_mark_time" => stats.total_mark_time/1e9,  # Convert to seconds
        "gc_sweep_time" => stats.total_sweep_time/1e9  # Convert to seconds
    )
end

function run(L::Int,p_ctrl::Float64,p_proj::Float64,seed::Int,ancilla::Int,maxdim::Int)
    println("Initializing...")
    init_stats = get_mem_stats()
    init_time = @elapsed begin
        ct_f=CT.CT_MPS(L=L,seed=seed,folded=true,store_op=false,store_vec=false,ancilla=ancilla,debug=false,xj=Set([1//3,2//3]),_maxdim=maxdim)
    end
    init_mem = get_mem_stats()
    println("Initialization took $(init_time) seconds")
    println("Memory allocated during init: $((init_mem["allocated_bytes"] - init_stats["allocated_bytes"])/1024/1024) MB")
    println("GC time during init: $(init_mem["gc_total_time"] - init_stats["gc_total_time"]) seconds")
    
    i=1
    T_max = ancilla ==0 ? 2*(ct_f.L^2) : div(ct_f.L^2,2)

    println("Starting main loop...")
    loop_stats = get_mem_stats()
    loop_time = @elapsed begin
        prog = Progress(T_max, 1, "Processing iterations...")
        for idx in 1:T_max
            if idx % 10 == 0  # Print memory stats every 10 iterations
                curr_stats = get_mem_stats()
                mem_allocated = (curr_stats["allocated_bytes"] - loop_stats["allocated_bytes"])/1024/1024
            end
            i=CT.random_control!(ct_f,i,p_ctrl,p_proj)
            next!(prog)
        end
    end
    final_stats = get_mem_stats()
    println("\nMain loop took $(loop_time) seconds")
    println("Memory allocated during loop: $((final_stats["allocated_bytes"] - loop_stats["allocated_bytes"])/1024/1024) MB")
    println("GC time during loop: $(final_stats["gc_total_time"] - loop_stats["gc_total_time"]) seconds")
    
    O=CT.order_parameter(ct_f)
    max_bond= CT.max_bond_dim(ct_f.mps)
    
    # Calculate total memory stats
    total_allocated = (final_stats["allocated_bytes"] - init_stats["allocated_bytes"])/1024/1024  # Convert to MB
    total_gc_time = final_stats["gc_total_time"] - init_stats["gc_total_time"]
    
    memory_stats = Dict(
        "total_allocated_mb" => total_allocated,
        "total_gc_time" => total_gc_time,
        "gc_pause_time" => final_stats["gc_pause_time"] - init_stats["gc_pause_time"],
        "gc_mark_time" => final_stats["gc_mark_time"] - init_stats["gc_mark_time"],
        "gc_sweep_time" => final_stats["gc_sweep_time"] - init_stats["gc_sweep_time"]
    )
    
    if ancilla ==0 
        EE=CT.von_Neumann_entropy(ct_f.mps,div(ct_f.L,2))
        return Dict(
            "O" => O, 
            "EE" => EE, 
            "max_bond" => max_bond, 
            "init_time" => init_time, 
            "loop_time" => loop_time, 
            "total_time" => (init_time + loop_time),
            "memory_stats" => memory_stats
        )
    else
        SA=CT.von_Neumann_entropy(ct_f.mps,1)
        return Dict(
            "O" => O, 
            "SA" => SA, 
            "max_bond" => max_bond, 
            "init_time" => init_time, 
            "loop_time" => loop_time, 
            "total_time" => (init_time + loop_time),
            "memory_stats" => memory_stats
        )
    end
end

function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--p_ctrl", "-c"
        arg_type = Float64
        default = 0.0
        help = "control map rate"
        "--p_proj", "-j"
        arg_type = Float64
        default = 0.0
        help = "projection rate"
        "--L", "-L"
        arg_type = Int
        default = 8
        help = "system size"
        "--seed", "-s"
        arg_type = Int
        default = 0
        help = "random seed"
        "--ancilla", "-a"
        arg_type = Int
        default = 0
        help = "number of ancilla"
        "--maxdim", "-m"
        arg_type = Int
        default = 10
        help = "set the maximal bond dim"
    end
    return parse_args(s)
end

function main()
    println("Uses threads: ",BLAS.get_num_threads())
    println("Uses backends: ",BLAS.get_config())
    args = parse_my_args()
    results = run(args["L"], args["p_ctrl"], args["p_proj"], args["seed"], args["ancilla"], args["maxdim"])

    filename = "MPS_(1-3,2-3)_L$(args["L"])_p_ctrl$(@sprintf("%.3f", args["p_ctrl"]))_p_proj$(@sprintf("%.3f", args["p_proj"]))_s$(args["seed"])_a$(args["ancilla"]).json"
    data_to_serialize = merge(results, Dict("args" => args))
    json_data = JSON.json(data_to_serialize)
    open(filename, "w") do f
        write(f, json_data)
    end
end

if isdefined(Main, :PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    main()
end




# windows: julia --sysimage ~/.julia/sysimages/sys_itensors.so run_CT_MPS_1-3.jl --p_ctrl 0.5 --p_proj 0.0 --L 8 --seed 0 --ancilla 0
# linux/mac: julia --sysimage ~/.julia/sysimages/sys_itensors.dylib run_CT_MPS_1-3.jl --p_ctrl 0.5 --p_proj 0.0 --L 8 --seed 0 --ancilla 0