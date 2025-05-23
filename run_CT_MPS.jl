using Pkg
using ITensors
using Random
using LinearAlgebra
using MKL
using JSON
using CT
using Printf
using ArgParse
using Serialization

function run(L::Int, p_ctrl::Float64, p_proj::Float64, seed::Int, ancilla::Int, maxdim::Int, store_trajectory::Bool=false)
    ct_f=CT.CT_MPS(L=L,seed=seed,folded=true,store_op=false,store_vec=false,ancilla=ancilla,debug=false,xj=Set([0]),_maxdim=maxdim)
    i=1
    T_max = ancilla ==0 ? 2*(ct_f.L^2) : div(ct_f.L^2,2)

    if store_trajectory
        O_list = zeros(T_max + 1)
        EE_list = zeros(T_max + 1)
        O_list[1] = CT.order_parameter(ct_f)
        EE_list[1] = ancilla == 0 ? CT.von_Neumann_entropy(ct_f.mps,div(ct_f.L,2)) : CT.von_Neumann_entropy(ct_f.mps,1)
    end

    for idx in 1:T_max
        println(idx)
        i=CT.random_control!(ct_f,i,p_ctrl,p_proj)
        if store_trajectory
            O_list[idx + 1] = CT.order_parameter(ct_f)
            EE_list[idx + 1] = ancilla == 0 ? CT.von_Neumann_entropy(ct_f.mps,div(ct_f.L,2)) : CT.von_Neumann_entropy(ct_f.mps,1)
        end
    end

    O = CT.order_parameter(ct_f)
    max_bond = CT.max_bond_dim(ct_f.mps)
    
    if store_trajectory
        return Dict("O" => O_list, "EE" => EE_list, "max_bond" => max_bond)
    else
        if ancilla == 0 
            EE = CT.von_Neumann_entropy(ct_f.mps,div(ct_f.L,2))
            return Dict("O" => O, "EE" => EE, "max_bond" => max_bond)
        else
            SA = CT.von_Neumann_entropy(ct_f.mps,1)
            return Dict("O" => O, "SA" => SA, "max_bond" => max_bond)
        end
    end
end

function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--p_ctrl", "-c"
        arg_type = Float64
        default = 0.0
        help = "measurement rate"
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
        "--ancilla", "-a"
        arg_type = Int
        default = 0
        help = "number of ancilla"
        "--maxdim", "-m"
        arg_type = Int
        default = 10
        help = "set the maximal bond dim"
        "--store_trajectory", "-t"
        action = :store_true
        help = "store the full time evolution trajectory"
    end
    return parse_args(s)
end

function main()
    println("Uses threads: ",BLAS.get_num_threads())
    println("Uses backends: ",BLAS.get_config())
    args = parse_my_args()
    results = run(args["L"], args["p_ctrl"], args["p_proj"], args["seed"], args["ancilla"], args["maxdim"], args["store_trajectory"])

    filename = if args["store_trajectory"]
        "MPS_(0,1)_L$(args["L"])_pctrl$(@sprintf("%.3f", args["p_ctrl"]))_pproj$(@sprintf("%.3f", args["p_proj"]))_s$(args["seed"])_T.json"
    else
        "MPS_(0,1)_L$(args["L"])_c$(@sprintf("%.3f", args["p_ctrl"]))_p$(@sprintf("%.3f", args["p_proj"]))_s$(args["seed"])_a$(args["ancilla"]).json"
    end
    
    data_to_serialize = merge(results, Dict("args" => args))
    json_data = JSON.json(data_to_serialize)
    open(filename, "w") do f
        write(f, json_data)
    end
end

if isdefined(Main, :PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    main()
end

# julia --sysimage ~/.julia/sysimages/sys_itensors.so run_CT_MPS.jl --p_ctrl 0.5 --p_proj 0. --L 8 --seed 0 --ancilla 0