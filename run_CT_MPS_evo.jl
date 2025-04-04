using ITensors
using Random
using LinearAlgebra
using MKL
using Pkg
using JSON
Pkg.activate("CT")
using CT
using Printf

using ArgParse
using Serialization
""" compute domain wall as a function of t"""

function random_int(L,lower_bound,upper_bound,seed=nothing)
    # lower_bound = 2^(L-1)
    # upper_bound = 2^L - 1
    if seed !== nothing
        rng = MersenneTwister(seed)
        return rand(rng, lower_bound:upper_bound)
    else
        return rand(lower_bound:upper_bound)
    end
end

function run_dw_t(L::Int,p_ctrl::Float64,p_proj::Float64,seed_C::Int,seed_m::Int,filename,args,)
    ct=CT.CT_MPS(L=L,seed=0,seed_C=seed_C,seed_m=seed_m,folded=true,store_op=true,store_vec=false,ancilla=0,xj=Set([0]),x0=1//2^L)
    print("x0: ", ct.x0)
    i=L
    tf=(ct.ancilla ==0) ? 2*ct.L^2 : div(ct.L^2,2)
    maxbond = CT.max_bond_dim(ct.mps)
    success = true
    for idx in 1:tf
        println(idx)
        try
            i=CT.random_control!(ct,i,p_ctrl,p_proj)
        catch e
            # Handle the error
            println("Caught an error: ", e)
            success = false
            break
        finally
            # Code that always executes, regardless of error
            maxbond = CT.max_bond_dim(ct.mps)
            success = idx == tf ? true : false
            results = Dict("maxbond"=>maxbond,"success"=>success,"idx"=>idx)
            open(filename, "w") do f
                data_to_serialize = merge(results, Dict("args" => args))
                json_data = JSON.json(data_to_serialize)
                write(f, json_data)
            end
        end
        
    end

    # O=CT.Z_bitstring(CT.bitstring_sample(ct))
    return Dict("maxbond"=>maxbond,"success"=>success)
end


function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--p_ctrl"
        arg_type = Float64
        default = 0.0
        help = "control rate"
        "--p_proj"
        arg_type = Float64
        default = 0.0
        help = "projection rate"
        "--L", "-L"
        arg_type = Int
        default = 8
        help = "system size"
        "--seed_C", "-C"
        arg_type = Int
        default = 0
        help = "random seed for circuit-- unitary and position of projection"
        "--seed_m", "-m"
        arg_type = Int
        default = 0
        help = "random seed for measurement outcome"
    end
    return parse_args(s)
end

function main()
    println("Uses threads: ",BLAS.get_num_threads())
    println("Uses backends: ",BLAS.get_config())
    args = parse_my_args()
    results = run_dw_t(args["L"], args["p_ctrl"], args["p_proj"], args["seed_C"],args["seed_m"])

    filename = "MPS_(0,1)_L$(args["L"])_pctrl$(@sprintf("%.3f", args["p_ctrl"]))_pproj$(@sprintf("%.3f", args["p_proj"]))_sC$(args["seed_C"])_sm$(args["seed_m"])_x01_evo.json"
    data_to_serialize = merge(results, Dict("args" => args))
    json_data = JSON.json(data_to_serialize)
    open(filename, "w") do f
        write(f, json_data)
    end
end


function main_interactive(L::Int,p_ctrl::Float64,p_proj::Float64,seed_C::Int,seed_m::Int)
    start_time = time()

    # println("Uses threads: ",BLAS.get_num_threads())
    # println("Uses backends: ",BLAS.get_config())
    # args = parse_my_args()
    args=Dict("L"=>L,"p_ctrl"=>p_ctrl,"p_proj"=>p_proj,"seed_C"=>seed_C,"seed_m"=>seed_m)
    filename = "MPS_(0,1)_L$(args["L"])_pctrl$(@sprintf("%.3f", args["p_ctrl"]))_pproj$(@sprintf("%.3f", args["p_proj"]))_sC$(args["seed_C"])_sm$(args["seed_m"])_x01_evo.json"
    
    results = run_dw_t(L, p_ctrl, p_proj, seed_C,seed_m,filename,args)
    
    elapsed_time = time() - start_time
    println("p_ctrl: ", args["p_ctrl"], " p_proj: ", p_proj, " L: ", L, " seed_C: ", seed_C, " seed_m: ", seed_m)
    println("Execution time: ", elapsed_time, " s")
end

if isdefined(Main, :PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    main()
end




# julia run_CT_MPS_evo.jl --p_ctrl 0.5 --p_proj 0.0 --L 8 --seed_C 0 --seed_m 0