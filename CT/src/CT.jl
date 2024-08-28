module CT

using ITensors
using Random
using LinearAlgebra
import TensorCrossInterpolation as TCI
using TCIITensorConversion

# using TimerOutputs
# const to = TimerOutput()
function __init__()
    ITensors.disable_warn_order()
end

mutable struct CT_MPS
    L::Int
    store_vec::Bool
    store_op::Bool
    store_prob::Bool
    seed::Union{Int,Nothing}
    seed_vec::Union{Int,Nothing}
    seed_C::Union{Int,Nothing} # random seed for "circuit", i.e., choice of unitary (Bernoulli) and position of measurement (control)
    seed_m::Union{Int,Nothing} # random seed for "measurement outcome".
    x0::Union{Rational{Int},Nothing}
    xj::Set
    _eps::Float64
    ancilla::Int
    folded::Bool
    rng::Random.AbstractRNG
    rng_vec::Random.AbstractRNG
    rng_C::Random.AbstractRNG
    rng_m::Random.AbstractRNG
    qubit_site::Vector{Index{Int64}}
    phy_ram::Vector{Int}
    ram_phy::Vector{Int}
    phy_list::Vector{Int}
    _maxdim::Int
    _cutoff::Float64
    mps::MPS
    vec_history::Vector{MPS}
    op_history::Vector{Vector{Any}}
    adder::Vector{Union{MPO,Nothing}}
    dw::Vector{Vector{Union{MPO,Nothing}}}
    debug::Bool
end

function CT_MPS(
    ; L::Int,
    store_vec::Bool=false,
    store_op::Bool=false,
    store_prob::Bool=false,
    seed::Union{Nothing,Int}=nothing,
    seed_vec::Union{Nothing,Int}=nothing,
    seed_C::Union{Nothing,Int}=nothing,
    seed_m::Union{Nothing,Int}=nothing,
    x0::Union{Rational{Int},Nothing}=nothing,
    xj::Set=Set([1 // 3, 2 // 3]),
    _eps::Float64=1e-10,
    ancilla::Int=0,
    folded::Bool=false,
    _maxdim::Int=10,
    _cutoff::Float64=1e-10,
    debug::Bool=false
    )
    rng = MersenneTwister(seed)
    rng_vec = seed_vec === nothing ? rng : MersenneTwister(seed_vec)
    rng_C = seed_C === nothing ? rng : MersenneTwister(seed_C)
    rng_m = seed_m === nothing ? rng : MersenneTwister(seed_m)
    qubit_site, ram_phy, phy_ram, phy_list = _initialize_basis(L,ancilla,folded)
    mps=_initialize_vector(L,ancilla,x0,folded,qubit_site,ram_phy,phy_ram,phy_list,rng_vec,_cutoff,_maxdim)
    adder=[adder_MPO(i1,xj,qubit_site,L,phy_ram,phy_list) for i1 in 1:L]
    dw=[[dw_MPO(i1,xj,qubit_site,L,phy_ram,phy_list,order) for i1 in 1:L] for order in 1:2]
    ct = CT_MPS(L, store_vec, store_op, store_prob, seed, seed_vec, seed_C, seed_m, x0, xj, _eps, ancilla, folded, rng, rng_vec, rng_C, rng_m, qubit_site, phy_ram, ram_phy, phy_list, _maxdim, _cutoff, mps, [],[],adder,dw,debug)
    return ct
end

function _initialize_basis(L,ancilla,folded)

    qubit_site = siteinds("Qubit", L+ancilla) # RAM site index
    # ram_phy[actual in ram] = physical 
    if ancilla ==0
        ram_phy = folded ? [i for pairs in zip(1:(L÷2), reverse((L÷2+1):L)) for i in pairs] : collect(1:L)
    elseif ancilla ==1
        ram_phy = folded ? vcat(L+1,[i for pairs in zip(1:(L÷2), reverse((L÷2+1):L)) for i in pairs]) : vcat(collect(1:L),L+1)
        # if folded : [L+1, 1,L, ..] else [1,2,..,L,L+1]
    elseif ancilla ==2
        error("Not implemented yet")
    end

    # phy_ram[physical] = actual in ram
    phy_ram = fill(0, L+ancilla)
    for (ram, phy) in enumerate(ram_phy)
        phy_ram[phy] = ram
    end

    # phy_list = Dict(0 => 1:L,1=>1:L)[ancilla]
    phy_list = collect(1:L)
    return qubit_site, ram_phy, phy_ram, phy_list
end

function _initialize_vector(L,ancilla,x0,folded,qubit_site,ram_phy,phy_ram,phy_list,rng_vec,_cutoff,_maxdim)
    if ancilla == 0
        if x0 !== nothing
            vec_int = dec2bin(x0, L)
            vec_int_pos = [string(s) for s in lpad(string(vec_int, base=2), L, "0")] # physical index
            return MPS(ComplexF64, qubit_site, [vec_int_pos[ram_phy[i]] for i in 1:L])
        else
            return randomMPS(rng_vec, qubit_site, linkdims=_maxdim)
        end
    elseif ancilla ==1
        if x0 !== nothing
            # This branch is only for debug purpose
            vec_int = dec2bin(x0, L)
            vec_int_pos = [string(s) for s in lpad(string(vec_int, base=2), L, "0")] # physical index
            push!(vec_int_pos,"0")
            return MPS(ComplexF64,qubit_site, [vec_int_pos[ram_phy[i]] for i in 1:(L+ancilla)])
        else
            # Create two orthogonal MPS and couple them to |0> and |1> ancilla, respectively
            # The physical qubit of ancilla is always the last one
            # The ram position, depending on the folded or not, is either the first (folded) or the last (unfolded)
            
            mps1 = randomMPS(rng_vec, ComplexF64, qubit_site[sort(phy_ram[phy_list])], linkdims=_maxdim)
            mps0 = randomMPS(rng_vec, ComplexF64, qubit_site[sort(phy_ram[phy_list])], linkdims=_maxdim)
            mps1 -= inner(mps0, mps1) .* mps0
            mps1 = mps1 ./ norm(mps1)
            anc0 = MPS(qubit_site[[(phy_ram[L+1])]], ["0"])
            anc1 = MPS(qubit_site[[(phy_ram[L+1])]], ["1"])
            if folded
                return add(attach_mps(anc0,mps0) , attach_mps(anc1,mps1),cutoff=_cutoff)/sqrt(2)
            else 
                return (attach_mps(mps0,anc0) + attach_mps(mps1,anc1))/sqrt(2)
            end
        end
    end
    # return MPS(vec, ct.qubit_site; cutoff=ct._eps, maxdim=ct._maxdim)
end

"""
apply the operator `op` to `mps`, the `op` should have indices of (i,j,k.. i',j',k')
the orthogonalization center is at i
index should be ram index
"""
function apply_op!(mps::MPS, op::ITensor, cutoff::Float64)
    i_list = [parse(Int, replace(string(tags(inds(op)[i])[length(tags(inds(op)[i]))]), "n=" => "")) for i in 1:div(length(op.tensor.inds), 2)]
    sort!(i_list)
    # println(i_list)
    orthogonalize!(mps, i_list[1])
    mps_ij = mps[i_list[1]]
    for idx in i_list[1]+1:i_list[end]
        mps_ij *= mps[idx]
    end
    mps_ij *= op 
    noprime!(mps_ij)
    
    if length(i_list) == 1
        mps[i_list[1]] = mps_ij
    else
        lefttags= (i_list[1]==1) ? nothing : tags(linkind(mps,i_list[1]-1))
        for idx in i_list[1]:i_list[end]-1
            inds1 = (idx ==1) ? [siteind(mps,1)] : [findindex(mps[idx-1],lefttags), findindex(mps[idx],"Site")]
            lefttags=tags(linkind(mps,idx))
            U, S, V = svd(mps_ij, inds1, cutoff=cutoff,lefttags=lefttags)
            mps[idx] = U
            mps_ij = S * V
        end
        mps[i_list[end]] = mps_ij
    end
    return
end
""" apply Bernoulli_map to physical site i
"""
function Bernoulli_map!(ct::CT_MPS, i::Int)
    S!(ct, i, ct.rng_C)
end

""" apply scrambler (Haar random unitary) to site (i,i+1) [physical index]
"""
function S!(ct::CT_MPS, i::Int, rng::Random.AbstractRNG; builtin=false)
    # U=ITensor(1.)
    # U *= randomUnitary(linkind(mps,i), linkind(mps,i+1))
    # mps[i] *= U
    # mps[i+1] *= U
    # return
    U_4_mat=U(4,rng)
    U_4 = reshape(U_4_mat, 2, 2, 2, 2)
    if ct.ancilla == 0 || ct.ancilla ==1
        ram_idx = ct.phy_ram[[ct.phy_list[i], ct.phy_list[(i)%(ct.L)+1]]]
        U_4_tensor = ITensor(U_4, ct.qubit_site[ram_idx[1]], ct.qubit_site[ram_idx[2]], ct.qubit_site[ram_idx[1]]', ct.qubit_site[ram_idx[2]]')
        # return U_4_tensor
        if builtin
            ct.mps=apply(U_4_tensor,ct.mps;cutoff=ct._cutoff)
        else
            apply_op!(ct.mps, U_4_tensor,ct._cutoff)
        end

        if ct.debug
            println("U $(U_4_mat) apply to $(i)")
        end
    elseif ct.ancilla ==2
        nothing
    end
end

""" apply control_map to physical sites i, i[1] happens to be the leading bit
"""
function control_map(ct::CT_MPS, n::Vector{Int}, i::Vector{Int})
    R!(ct, n, i)
    if ct.xj == Set([1 // 3, 2 // 3])
        ct.mps=apply(ct.adder[i[1]],ct.mps;cutoff=ct._cutoff)
        normalize!(ct.mps)
        truncate!(ct.mps, cutoff=ct._cutoff)
    end
end

"""Rest: apply Projection to physical site list i with outcome list n, then reset to 1
"""
function R!(ct::CT_MPS, n::Vector{Int}, i::Vector{Int})
    P!(ct, n, i)
    if ct.xj in Set([Set([1 // 3, 2 // 3]),Set([0])])
        if n[1]==1
            X!(ct,i[1])
            # print("X")
        end
    elseif ct.xj in [Set([1 // 3, -1 // 3])]
        if n[1]^n[2]==0
            print("i[end] needs to be checked")
            X!(ct,i[end])
        end
    end
end
"""generate projector for physical site (will convert to RAM site)
"""
function projector(ct::CT_MPS,n::Vector{Int}, i::Vector{Int})
    if ct.debug
        println("Get projector $(n) at Phy $(i)")
    end
    proj_op=emptyITensor(ct.qubit_site[ct.phy_ram[ct.phy_list[i]]],ct.qubit_site[ct.phy_ram[ct.phy_list[i]]]')
    idx=n.+1
    proj_op[ idx...,idx... ]=1+0im
    return proj_op
end
""" perform projection for physical site
"""
function P!(ct::CT_MPS, n::Vector{Int}, i::Vector{Int})
    @assert length(n) == length(i) "length of n $(n) is not equal to length of i $(i)"
    if ct.debug
        println("Projecting $(n) at Phy $(i)")
    end
    proj_op= projector(ct,n, i)
    apply_op!(ct.mps, proj_op, ct._cutoff)
    if ct.debug
        println("norm is $(norm(ct.mps))")
    end
    normalize!(ct.mps)
    truncate!(ct.mps, cutoff=ct._cutoff)
end
""" perform Sigma_X for physical site
"""
function X!(ct::CT_MPS, i::Int)
    X_op=ITensor([0 1+0im; 1+0im 0],ct.qubit_site[ct.phy_ram[ct.phy_list[i]]],ct.qubit_site[ct.phy_ram[ct.phy_list[i]]]')
    apply_op!(ct.mps, X_op, ct._cutoff)
    # normalize!(ct.mps)
end


"""
compute the entanglement entropy of 1...i, i+1...L (ram sites)
"""
function von_Neumann_entropy(mps::MPS, i::Int)
    mps_ = orthogonalize(mps, i)
    # _,S = svd(mps_[i], (linkind(mps_, i-1), siteind(mps_,i)))
    _, S = svd(mps_[i], (linkind(mps_, i),))
    SvN = 0.0
    # print(S)
    for n = 1:dim(S, 1)
        p = S[n, n]^2
        SvN -= p * log(p)
    end
    return SvN
end

function mps_to_tensor(mps; array=false, vector=false, column_first=true)
    # psi = mps[1]
    # for i = 2:length(mps)
    #     psi = psi * mps[i]
    # end

    psi=prod(mps)
    # Convert the resulting tensor to a dense array
    if array
        psi = array(psi)
        if vector
            if column_first
                return vec(psi)
            else
                return vec(permutedims(psi, reverse(1:ndims(psi))))
            end
        else
            return psi
        end
    else
        return psi
    end
end

# function encoding!(ct::CT_MPS)
#     i=1
#     for idx in 1:2*ct.L^2
#         Bernoulli_map!(ct, i)
#         # println(i)
#         i=mod(((i+1) - 1),ct.L+ct.ancilla )+ 1
#     end
# end
"""randomly apply control or Bernoulli map to physical site i (the left leg of op)
"""
function random_control!(ct::CT_MPS, i::Int, p_ctrl, p_proj)
    
    op_l=[]
    p_0=-1.  # -1 for not applicable because of Bernoulli map
    if rand(ct.rng_C) < p_ctrl
        # control map
        if ct.xj in Set([Set([1 // 3, 2 // 3]),Set([0])])
            p_0= inner_prob(ct, [0], [i])
            if ct.debug
                println("Born prob for measuring 0 at phy site $i is $p_0")
            end
            n =  rand(ct.rng_m) < p_0 ?  0 : 1
            control_map(ct, [n], [i])
            push!(op_l,Dict("Type"=>"Control","Site"=>[i],"Outcome"=>[n]))
        elseif ct.xj in [Set([1 // 3, -1 // 3])]
            # p_00= ...
            # p_01= ...
            # p_10= ...
            # p_11= ...
            # n = ...
            # control_map(ct, n, [i, i+1])
            nothing
        end
        if ct.debug
            print("Control with $(i)")
        end
        i=mod(((i-1) - 1),(ct.L)) + 1
        if ct.debug
            println("=> Next i is $(i)")
        end
    else
        # chaotic map
        Bernoulli_map!(ct, i)
        push!(op_l,Dict("Type"=>"Bernoulli","Site"=>[i,((i+1) - 1)%(ct.L) + 1],"Outcome"=>nothing))
        i=mod(((i+1) - 1),(ct.L) )+ 1
    end

    if op_l[end]["Type"] == "Bernoulli"
        # projection
        for pos in [i-1,i]
            if rand(ct.rng_C) < p_proj
                pos=mod((pos-1),ct.L)+1
                p2=inner_prob(ct, [0], [pos])
                n= rand(ct.rng_m) < p2 ? 0 : 1
                P!(ct,[n],[pos])
                push!(op_l,Dict("Type"=>"Projection","Site"=>[pos],"Outcome"=>[n]))
            end
        end
    end
    update_history(ct,op_l,p_0)
    
    return i
end

function update_history(ct::CT_MPS,op::Vector{Any},p_0::Float64)
    if ct.store_vec
        push!(ct.vec_history,copy(ct.mps)) 
    end
    if ct.store_op
        push!(ct.op_history,op)
    end
    if ct.store_prob
        push!(ct.op_history,p_0)
    end
end
    

""" compute the Born probability through the inner product at physical site list i (will convert to RAM site internally)
"""
function inner_prob(ct::CT_MPS, n::Vector{Int}, i::Vector{Int})
    @assert length(n) == length(i) "length of n $(n) is not equal to length of i $(i)"
    ram_idx=ct.phy_ram[ct.phy_list[i]]
    if length(i)==1
        if ct.debug
            println("Get projector for inner_prob $(n) at Phy $(i) at RAM $(ram_idx)")
        end
        proj_op= array(projector(ct,n, i))
        return only(expect(ct.mps, proj_op,sites=ram_idx))
    else
        # proj_op= projector(ct,n, i,tensor=false)
        # os=OpSum()
        error("Not implemented yet")
    end

end

function order_parameter(ct::CT_MPS)
    if ct.xj in Set([Set([1 // 3, 2 // 3])])
        O = ZZ(ct)
        return (O)
    elseif ct.xj in [Set([0])]
        O = Z(ct)
        return (O)
    end
end

function Zi(ct::CT_MPS)
    sZ = expect(ct.mps, "Sz")[ct.phy_ram][ct.phy_list]
    return sZ*2
end

"""note that here ZiZj is not put to the correct position because it would be taken a sum in the end anyway""" 
function ZiZj(ct::CT_MPS)
    sZZ=correlation_matrix(ct.mps, "Sz", "Sz")
    return sZZ*4
end

function Z(ct::CT_MPS)
    sZ = Zi(ct)
    if ct.debug
        println("Z are $(2 * sZ)")
    end
    return real(sum(sZ)) / ct.L
end

function Z_sq(ct::CT_MPS)
    sZZ = ZiZj(ct)
    return real(sum(sZZ)) / ct.L^2
end

function ZZ(ct::CT_MPS)
    os = OpSum()
    for i in ct.phy_list[1:end-1]
        os += "Sz", ct.phy_ram[i], "Sz", ct.phy_ram[i+1]
    end
    os += "Sz", ct.phy_ram[ct.phy_list[end]], "Sz", ct.phy_ram[ct.phy_list[1]]
    zz = MPO(os, ct.qubit_site)
    return real(-inner(ct.mps', zz, ct.mps)) * 4 / ct.L
end

function mps_element(mps::MPS, x::String)
    @assert length(x) == length(mps)
    x = [string(s) for s in x]
    V = ITensor(1.0)
    for i = 1:length(mps)
        V *= (mps[i] * state(siteind(mps, i), x[i]))
    end
    return scalar(V)
end

function display_mps_element(ct::CT_MPS)
    println(rpad("RAM",ct.L+ct.ancilla), "=>", "Physical")
    vec=zeros(Complex{Float64},2^(ct.L+ct.ancilla))
    for i in 1:2^(ct.L+ct.ancilla)
        bitstring=lpad(string(i-1,base=2),ct.L+ct.ancilla,"0")
        matel=CT.mps_element(ct.mps,bitstring)
        println(bitstring, "=>", bitstring[ct.phy_ram],": ",matel)
        vec[parse(Int,bitstring[ct.phy_ram];base=2)+1]=matel
    end
    return vec
end

function dec2bin(x::Real, L::Int)
    @assert 0 <= x < 1 "$x is not in [0,1)"

    return BigInt(floor(x * (BigInt(1) << L)))
end

"""create Haar random unitary
"""
function U(n, rng::Random.AbstractRNG=MersenneTwister(nothing))
    z = randn(rng, 4, 4) + randn(rng, 4, 4) * im
    Q, R = qr(z)
    r_diag = diag(R)
    Lambda = Diagonal(r_diag ./ abs.(r_diag))
    Q *= Lambda
    return Q
end
"""attach mps2 to mps1 at the last site. The returned mps is |mps1>|mps2>
The physical index of mps2 should be given 
"""
function attach_mps(mps1::MPS,mps2::MPS,cutoff::Float64=1e-10;replace_link=true)
    n1=length(mps1)
    n2=length(mps2)
    if replace_link
        advance_link_tags!(mps2,n1)
    end
    orthogonalize!(mps1,n1)
    mps_=mps1[n1]*mps2[1]
    lefttags=replace(string(tags(findinds(inds(mps_),"Link")[1])),r"l=[\d]+"=>"l=$(n1)",'"'=>"")
    U,S,V=svd(mps_,inds(mps1[n1]),cutoff=cutoff,lefttags=lefttags)
    mps_new=MPS(vcat(siteinds(mps1),siteinds(mps2)))
    for i in 1:(n1+n2)
        if i<n1
            mps_new[i]=mps1[i]
        elseif i==n1
            mps_new[i]=U
        elseif i==n1+1
            mps_new[i]=S*V
        else
            mps_new[i]=mps2[i-n1]
        end
    end
    return mps_new
end
function advance_link_tags!(mps::MPS,n::Int64)
    for i in 1:length(mps)
        links=findinds(mps[i],"Link")
        for link in links
            cur_i=parse(Int,match(r"l=(\d+)",string(tags(link)[length(tags(link))]))[1])
            new_link=replacetags(link,"l=$(cur_i)","l=$(n+cur_i)")
            mps[i]*=delta(link,new_link)
        end
    end
end

function max_bond_dim(mps::MPS)
    max_dim = 0
    for i in 1:length(mps)-1
        dim = commonind(mps[i], mps[i+1])
        max_dim = max(max_dim, space(dim))
    end
    return max_dim
end
"""add one, i1 is the leading (physical) qubit"""
function add1(i1::Int,L::Int,phy_ram::Vector{Int},phy_list::Vector{Int})
    A1=OpSum()
    A1+=tuple(([((i==i1) ? "X" : "S+",phy_ram[phy_list[i]]) for i in 1:L]...)...)

    for j in 2:L
        i_list = phy_list[mod.(collect(j:L).+(i1-1).-1,L).+1]
        A1+=tuple(([((i==i_list[1]) ? "S-" : "S+",phy_ram[i]) for i in i_list]...)...)
    end
    return A1
end
"""return domain wall position for FM fixed point
the dw = L P_{1}^1 + (L-1) P_{1}^0 P_{2}^1 + .. +  P_{1}^{0}P_{2}^{0}...P_{L}^{1}
To test: 
product state vs general state [x]
folded vs unfolded [x]
i1=1 vs i1 >1 [x]
---
order : the order of moment
"""
function dw_FM(i1::Int,L::Int,phy_ram::Vector{Int},phy_list::Vector{Int},order::Int)
    dw=OpSum()
    for j in 1:L
        i_list = phy_list[mod.(collect(1:j).+(i1-1).-1,L).+1]
        dw+=tuple((L-j+1)^order,([((i==i_list[end]) ? "Proj1" : "Proj0",phy_ram[i]) for i in i_list]...)...)
    end
    return dw
end

function dw_MPO(i1::Int,xj::Set,qubit_site::Vector{Index{Int64}},L::Int,phy_ram::Vector{Int},phy_list::Vector{Int},order::Int)
    if xj == Set([1 // 3, 2 // 3])
        return nothing
    elseif xj == Set([0])
        return MPO(dw_FM(i1,L,phy_ram,phy_list,order),qubit_site)
    end
end

function dw(ct::CT_MPS,i1::Int)
    if ct.xj == Set([1 // 3, 2 // 3])
        return nothing
    elseif ct.xj == Set([0])
        dw1=real(inner(ct.mps',ct.dw[1][i1],ct.mps))
        dw2=real(inner(ct.mps',ct.dw[2][i1],ct.mps))
        return dw1,dw2
    end
end

"""compute mpo^n
convert n to binary rep, and then square it each step, if the binary bit is 1, add the current mpo to the sum
the space overhead is O(2) forone n,
time complexity is O(L): (L-2)+(L/2)'s multiplication for 1/6, and (L-1)+ (L/2)'s multiplication for 1/3
"""
function power_mpo(mpo::MPO,n_list::Vector{Int})
    # mpo_sum=1
    L=length(mpo)
    n_list=mod.(n_list,2^L)
    mpo_sum=fill(copy(mpo),length(n_list))
    mpo_sum_visted=fill(false,L)
    # print(n_list)
    while sum(n_list)>0
        for (n_idx,n) in enumerate(n_list)
            if n%2==1
                # println(n_idx,':',n," in ",n_list)
                if !mpo_sum_visted[n_idx]
                    mpo_sum_visted[n_idx]=true
                    mpo_sum[n_idx]=mpo
                else
                    mpo_sum[n_idx]=apply(mpo_sum[n_idx],mpo)
                end
            end
        end
        mpo=apply(mpo,mpo)
        n_list.>>=1
    end
    return mpo_sum
end

"""return the MPO for {1/6,1/3}, here site indices are all physical indices"""
function adder_MPO(i1::Int,xj::Set,qubit_site::Vector{Index{Int64}},L::Int,phy_ram::Vector{Int},phy_list::Vector{Int})
    if xj == Set([1 // 3, 2 // 3])
        add1_mpo=MPO(add1(i1,L,phy_ram,phy_list),qubit_site)
        # print(add1_mpo)
        add1_6,add1_3=power_mpo(add1_mpo,[div(2^L,6)+1,div(2^L,3)])
        i2=phy_list[mod(i1,L)+1]    # 2
        add_condition=apply(add1_6,P_MPO([phy_ram[i2]],[0],qubit_site)) + apply(add1_3,P_MPO([phy_ram[i2]],[1],qubit_site))
        iLm2=phy_list[mod(i1+L-4,L)+1]  # L-2
        iLm1=phy_list[mod(i1+L-3,L)+1 ]   # L-1
        iL=phy_list[mod(i1+L-2,L)+1] # L i1+(L-1) -> L, (x-1)%L+1
        P2=(P_MPO([phy_ram[i1],phy_ram[iLm2],phy_ram[iL]],[1,0,1],qubit_site)+P_MPO([phy_ram[i1],phy_ram[iLm2],phy_ram[iL]],[0,1,0],qubit_site))    # fix spurs
        XI=XI_MPO([phy_ram[iLm1]],qubit_site)

        fix_spurs = apply(XI,P2) + I_MPO([phy_ram[iLm1]],qubit_site)
        return apply(fix_spurs,add_condition)
    else
        return nothing
    end
end
"""return the MPO for a joint projection at sites `pos_list` with outcome `n_list`"""
function P_MPO(pos_list::Vector{Int},n_list::Vector{Int},qubit_site::Vector{Index{Int64}})
    P_sum=OpSum()
    P_sum+=tuple(([("Proj$n",pos) for (pos,n) in zip(pos_list,n_list)]...)...)
    return MPO(P_sum,qubit_site)
end

"""return the oerpator of X-I at sites `pos_list`"""
function XI_MPO(pos_list::Vector{Int},qubit_site::Vector{Index{Int64}})
    X_sum=OpSum()
    X_sum+=tuple(([("X",pos) for pos in pos_list]...)...)
    X_sum+=tuple(([(-1,"I",pos) for pos in pos_list]...)...)
    return MPO(X_sum,qubit_site)
end

function I_MPO(pos_list::Vector{Int},qubit_site::Vector{Index{Int64}})
    I_sum=OpSum()
    I_sum+=tuple(([("I",pos) for pos in pos_list]...)...)
    return MPO(I_sum,qubit_site)
end

"""ket_index, bra_index, -1 for no projection, 0 for 0, 1 for 1, using physical index, i1 is the last qubit"""
function get_reduced_DM(dm::MPO,ct::CT_MPS,ket_index::Vector{Int},bra_index::Vector{Int},i1::Int)
    rdm=MPO(ct.qubit_site)
    for phy_idx in 1:ct.L
        ram_idx= ct.phy_ram[mod(ct.phy_list[phy_idx] + i1-1  ,ct.L)+1]
        # ram_idx= ct.phy_ram[ct.phy_list[phy_idx]]
        # println(mod(ct.phy_list[phy_idx] + i1 -1  ,ct.L)+1)
        ket_leg = ct.qubit_site[ram_idx]
        bra_leg = ket_leg'
        if ket_index[phy_idx] !=  -1
            rdm[ram_idx]=dm[ram_idx]* state(ket_leg,ket_index[phy_idx]+1)
        else
            rdm[ram_idx]=dm[ram_idx]
        end
        if bra_index[phy_idx] !=  -1
            rdm[ram_idx]=rdm[ram_idx]* state(bra_leg,bra_index[phy_idx]+1)
        end
    end
    truncate!(rdm,cutoff=ct._cutoff)
    return rdm
end

function get_DM(mps::MPS)
    return outer(mps',mps)
end

function l1_coherence_2(rho::MPO,ct::CT_MPS,k1::Int,k2::Int,i1::Int)
    L=length(rho)
    if k1==0
        ket_index=fill(0,L)
    else
        ket_index=vcat(fill(0, L - k1), [1], fill(-1, k1 - 1))
    end

    if k2==0
        bra_index = fill(0,L)
    else
        bra_index = vcat(fill(0, L - k2), [1], fill(-1, k2 - 1))
    end
    rdm = get_reduced_DM(rho,ct, ket_index, bra_index,i1)
    if k1 ==k2
        sum_=sum_of_norm(rdm,false)
        tr_rho = trace(rdm)
        intra_sum= (sum(length,siteinds(rdm)) < 44) ? sum_-tr_rho : sum_
        return  intra_sum, tr_rho
    else
        return sum_of_norm(rdm,true)
    end
end

function l1_coherence(rho::MPO,ct::CT_MPS,k1::Int,k2::Int,i1::Int)
    L=length(rho)
    if k1==0
        ket_index=fill(0,L)
    else
        ket_index=vcat(fill(0, L - k1), [1], fill(-1, k1 - 1))
    end

    if k2==0
        bra_index = fill(0,L)
    else
        bra_index = vcat(fill(0, L - k2), [1], fill(-1, k2 - 1))
    end
    rdm = get_reduced_DM(rho,ct, ket_index, bra_index,i1)
    if k1 ==k2
        sum_=sum_of_norm_MPS(rdm)
        tr_rho = trace(rdm)
        intra_sum= sum_-tr_rho
        return  intra_sum, tr_rho
    else
        return sum_of_norm_MPS(rdm)
    end
end

function sum_of_norm(rdm::MPO,inter)
    len_rdm= sum(length,siteinds(rdm)) 
    if len_rdm < 44
        return sum_of_norm_tensor(rdm)
    else
        return sum_of_norm_sample(rdm,inter)
    end
end

"""convert to tensor and sum"""
function sum_of_norm_tensor(rdm::MPO)
    return sum(abs.(prod(rdm)))
end

"""inner product with all one MPS (note that this needs to take a positive matrix, where every elemnt has bee taken norm already)"""
function sum_of_norm_MPS(rdm::MPO)
    L = length(rdm)
    site_idx=siteinds(rdm)
    V=ITensor(1.)
    for i in 1:L
        V*=rdm[i]
        for ind in site_idx[i]
            V*=ITensor([1.,1.],ind)
        end
    end
    return scalar(V)
end


function trace(rdm::MPO)
    for i in 1:length(rdm)
        if hastags(rdm[i],"Site")
            rdm[i] = tr(rdm[i])
        end
    end
    return abs(scalar(prod(rdm)))
end

function sum_of_norm_loop(rdm::MPO)
    sum_=0
    L= sum(length,siteinds(rdm)) 
    site_idx=siteinds(rdm)
        
    for basis in 1:2^L
        x=collect(lpad(string(basis-1,base=2),L,"0"))
        V = ITensor(1.0)
        for i in 1:length(rdm)
            V *= rdm[i] 
            for ind in site_idx[i]
                V *= state(ind,pop!(x)-'0'+1 )
            end
        end
        sum_ += abs(scalar(V))
    end
    return sum_
end

""" instead of summing all, using O(L) samples"""
function sum_of_norm_sample(rdm::MPO,inter::Bool)
    sum_=0
    L= sum(length,siteinds(rdm)) 
    site_idx=siteinds(rdm)
    basis_list = 2:L*2
    for basis in basis_list
        x=collect(lpad(string(basis-1,base=2),L,"0"))
        V = ITensor(1.0)
        for i in 1:length(rdm)
            V *= rdm[i] 
            for ind in site_idx[i]
                V *= state(ind,pop!(x)-'0'+1 )
            end
        end
        sum_ += abs(scalar(V)) 
    end
    if inter
        return sum_* 2^L/length(basis_list)
    else
        return sum_* (2^L-2^(L/2))/length(basis_list)
    end
end

function abs_mps(mps::MPS;tolerance::Real=1e-8, maxbonddim::Int=30)
    L=length(mps)
    sites=siteinds(mps)
    mps_tci=TCI.TensorTrain(mps)
    localdims = fill(2, L)
    mps_abs_func(v) = abs(mps_tci(v))
    mps_abs_func_cache = TCI.CachedFunction{Float64}(mps_abs_func,localdims)
    initialpivots=[fill(1,L),[1,2,fill(1,L-2)...]] # this gives the max value of wave function
    tci, _, _ = TCI.crossinterpolate2(Float64, mps_abs_func_cache, localdims, initialpivots; tolerance=tolerance, maxbonddim=maxbonddim,)
    return MPS(tci,sites=sites)
end

function get_coherence_matrix(ct::CT_MPS,i1::Int;tolerance::Real=1e-8, maxbonddim::Int=30)
    mps_abs= abs_mps(ct.mps;tolerance=tolerance,maxbonddim=maxbonddim)
    rho = get_DM(mps_abs)
    L=length(rho)
    coherence_matrix=zeros(L+1,L+1)
    fdw=zeros(L+1)
    for i in 0:L
        for j in 0:L
            if i == j
                # @timeit to "same" begin
                coherence_matrix[i+1,j+1], fdw[i+1] = l1_coherence(rho,ct,i,j,i1)
                # end
            else
                # @timeit to "different" begin
                coherence_matrix[i+1,j+1] = l1_coherence(rho,ct,i,j,i1)
                coherence_matrix[j+1,i+1] = coherence_matrix[i+1,j+1]
                # end
            end
        end
    end
    return coherence_matrix, fdw
end

"""return the element-wise product of two MPS mps1, and mp2"""
function elementwise_product(mps1::MPS, mps2::MPS;cutoff::Float64=1e-10,maxdim::Int=25,method::String="densitymatrix")
    site_idx=siteinds(mps1)
    mpo1=MPO(site_idx)
    for i in 1:length(mps1)
        mpo1[i]=mps1[i]* delta(site_idx[i],prime(site_idx[i],1),prime(site_idx[i],2))
    end
    mps_prod=apply(mpo1,prime(mps2;tags="Site"),method=method,cutoff=cutoff,maxdim=maxdim)
    # truncate!(mps_prod, cutoff=cutoff)
    return noprime!(mps_prod)
end

function all_one_mps(sites::Vector{Index{Int64}})
    mps=MPS(sites)
    for i in 1:length(mps)
        mps[i]=ITensor([1,1],sites[i])
    end
    return mps
end


"""return the element-wise sqrt of mps"""
function sqrt_mps(mps::MPS;eps::Float64=1e-5,maxiter::Int=100,cutoff::Float64=1e-10,maxdim::Int=20)
    L=siteinds(mps)
    allone=CT.all_one_mps(siteinds(mps))
    a=mps
    c=mps-allone
    diff_a = 100
    while diff_a > eps && maxiter > 0
        println(maxiter)
        a_new = -(a,elementwise_product(a,c,cutoff=cutoff,maxdim=maxdim)/2,cutoff=cutoff,maxdim=maxdim)
        truncate!(a_new)
        c = elementwise_product(elementwise_product(c,c, cutoff=cutoff,maxdim=maxdim),-(c,3*allone,cutoff=cutoff,maxdim=maxdim), cutoff=cutoff,maxdim=maxdim)/4
        truncate!(c)
        # diff_a=(1-inner(c,c))
        # diff_a=10
        diff_a_mps=a-a_new
        diff_a = abs(inner(conj(diff_a_mps),diff_a_mps))
        println(diff_a)
        a=a_new
        # c=c_new
        maxiter -= 1
    end
    return a,c
end


"""return the element-wise sqrt of mps, using 1/sqrt(S), potential caveat: divergence if S is too large"""
function sqrt_mps_inverse(mps::MPS,eps::Float64=1e-5,maxiter::Int=100)
    allone=CT.all_one_mps(siteinds(mps))
    x=allone
    diff_x = 100
    while diff_x > eps && maxiter > 0
        println(maxiter)
        x_new=elementwise_product(allone*3/2-elementwise_product(elementwise_product(x,x),mps)/2,x)
        diff_x_mps=x-x_new
        diff_x = inner(diff_x_mps,diff_x_mps)
        x=x_new
        maxiter -= 1
    end
    return elementwise_product(x,mps)
end

    

function test_profiler()
    ct=CT_MPS(L=4,seed=1,folded=true,store_op=false,store_vec=false,ancilla=0,_maxdim=6,xj=Set([0]),)
    dm = get_DM(ct.mps)
    coherence_matrix, fdw=get_coherence_matrix(dm,ct)
    # println("Timer outputs:")
    # show(to)
end

# function all_one_mps()


greet() = print("Hello World! How are 11?")


end # module CT
