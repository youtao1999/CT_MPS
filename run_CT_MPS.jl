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

function run(L::Int,p::Float64,seed::Int,ancilla::Int,maxdim::Int)
    
    ct_f=CT.CT_MPS(L=L,seed=seed,folded=true,store_op=false,store_vec=false,ancilla=ancilla,debug=false,xj=Set([0]),_maxdim=maxdim)
    i=1
    T_max = ancilla ==0 ? 2*(ct_f.L^2) : div(ct_f.L^2,2)

    for idx in 1:T_max
        println(idx)
        i=CT.random_control!(ct_f,i,p)
    end
    O=CT.order_parameter(ct_f)
    max_bond= CT.max_bond_dim(ct_f.mps)
    if ancilla ==0 
        EE=CT.von_Neumann_entropy(ct_f.mps,div(ct_f.L,2))
        return Dict("O" => O, "EE" => EE, "max_bond" => max_bond)
    else
        SA=CT.von_Neumann_entropy(ct_f.mps,1)
        return Dict("O" => O, "SA" => SA, "max_bond" => max_bond)
    end
end

function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--p", "-p"
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
    end
    return parse_args(s)
end

function main()
    println("Uses threads: ",BLAS.get_num_threads())
    println("Uses backends: ",BLAS.get_config())
    args = parse_my_args()
    results = run(args["L"], args["p"], args["seed"],args["ancilla"],args["maxdim"])

    # filename = "MPS_(0,1)_L$(args["L"])_p$(round(args["p"], digits=2))_s$(args["seed"]).jls"
    # open(filename, "w") do f
    #     serialize(f, Dict("O" => results[1], "EE" => results[2], "max_bond" => results[3],"args" => args))
    # end
    filename = "MPS_(0,1)_L$(args["L"])_p$(@sprintf("%.3f", args["p"]))_s$(args["seed"])_a$(args["ancilla"]).json"
    data_to_serialize = merge(results, Dict("args" => args))
    # Dict("O" => results[1], "EE" => results[2], "max_bond" => results[3], "args" => args)
    json_data = JSON.json(data_to_serialize)
    open(filename, "w") do f
        write(f, json_data)
    end
end

if isdefined(Main, :PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    main()
end




# julia --sysimage ~/.julia/sysimages/sys_itensors.so run_CT_MPS.jl --p 1 --L 8 --seed 0 --ancilla 0