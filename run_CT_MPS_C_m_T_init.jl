using Base.Threads

# p_ctrl_list = [0.4,0.5,0.6]
p_ctrl_list = [0.4,0.5,0.5]
p_proj=0.
L_list= [10]
# sC_list=range(200,999)
sC_list=range(200,500)
sm_list=range(0,499)
# sm_list=range(0,0)


include("run_CT_MPS_C_m_T.jl")
for p_ctrl in p_ctrl_list
    for L in L_list
        for sC in sC_list
            # Threads.@threads for sm in sm_list
            for sm in sm_list
                # println("p_ctrl: ", p_ctrl, " L: ", L, " sC: ", sC, " sm: ", sm)
                # ARGS = ["--p_ctrl", string(p_ctrl), "--p_proj", string(p_proj), "--L", string(L), "--seed_C", string(sC), "--seed_M", string(sm)]
                
                main_interactive(L,p_ctrl,p_proj,sC,sm)
            end
        end
    end
end


# julia run_CT_MPS_C_m_T.jl --p_ctrl 0.5 --p_proj 0.0 --L 10 --seed_C 200 --seed_m 1