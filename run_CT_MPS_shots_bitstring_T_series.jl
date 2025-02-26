include("run_CT_MPS_shots_bitstring_T.jl")

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
    for (p_ctrl,p_proj,L,seed_C,seed_m) in Iterators.partition(args, 5)
        main_interactive(parse(Int,L),parse(Float64,p_ctrl),parse(Float64,p_proj),parse(Int,seed_C),parse(Int,seed_m))
    end
end

if isdefined(Main, :PROGRAM_FILE) && abspath(PROGRAM_FILE) == @__FILE__
    main()
end


# julia run_CT_MPS_shots_bitstring_T_series.jl --params "0.5,0.0,8,0,0,0.5,0.0,8,0,1"