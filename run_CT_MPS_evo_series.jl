include("run_CT_MPS_evo.jl")

function parse_my_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--params"
        arg_type = String
        help = "Comma-separated list of values"
    end
    return parse_args(s)
end

function main()
    args = parse_my_args()
    args = split(args["params"],",")
    for (p_ctrl,p_proj,L,seed_C,seed_m,maxdim,cutoff) in Iterators.partition(args, 7)
        println("p_ctrl: ", p_ctrl, " p_proj: ", p_proj, " L: ", L, " seed_C: ", seed_C, " seed_m: ", seed_m, " maxdim: ", maxdim, " cutoff: ", cutoff)
        main_interactive(parse(Int,L),parse(Float64,p_ctrl),parse(Float64,p_proj),parse(Int,seed_C),parse(Int,seed_m),parse(Int,maxdim),parse(Float64,cutoff))
    end
end

if isdefined(Main, :PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    main()
end


# julia run_CT_MPS_evo_series.jl --params "0.1,0.0,8,0,0,10000,1e-10,0.1,0.0,8,0,1,10000,1e-10"
# julia run_CT_MPS_evo_series.jl --params "0.1,0.0,30,0,0,0.1,0.0,30,0,1"