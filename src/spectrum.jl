# AffineMapDiffEq.jl Spectrum and Eigen Analysis

"""
    spectrum(dp::dynamic_problemSampled; p=dp.Problem.p)

Calculates the Floquet multipliers (spectrum) of the mapping.
"""
function spectrum(dp::dynamic_problemSampled; p=dp.Problem.p)
    Nstep = size(dp.StateSmaplingTime, 1)
    s_start = randsimilar(dp.Problem.u0, Nstep)
    
    mus = getindex(schursolve(s -> LinMap(dp, s; p=p)[1], s_start, dp.Krylov_arg...), [3, 2, 1])
    return mus[1], mus[2]
end

"""
    spectralradius(dp::dynamic_problemSampled; p=dp.Problem.p)

Calculates the maximum absolute Floquet multiplier.
"""
function spectralradius(dp::dynamic_problemSampled; p=dp.Problem.p)
    if dp.zerofixpont
        return Float64(maximum(abs.(spectrum(dp; p=p)[1])))::Float64
    else
        return Float64(maximum(abs.(affine(dp; p=p)[1][1])))::Float64
    end
end

"""
    issi_eigen(dp::dynamic_problemSampled; p=dp.Problem.p)

Internal Subspace Iteration (ISSI) based eigen analysis.
"""
function issi_eigen(dp::dynamic_problemSampled; p=dp.Problem.p)
    Nstep = size(dp.StateSmaplingTime, 1)
    s0 = [zero(dp.Problem.u0) for _ in 1:Nstep]
    v0 = LinMap(dp, s0; p=p)[1]
    
    S = [randsimilar(dp.Problem.u0, Nstep) for _ in 1:dp.eigN]
    H = zeros(Float64, dp.eigN, dp.eigN)
    
    for _ in 1:12
        V = [LinMap(dp, Si + s0; p=p)[1] - v0 for Si in S]
        StS = [S[i]' * S[j] for i in 1:length(S), j in 1:length(S)]
        StV = [S[i]' * V[j] for i in 1:length(S), j in 1:length(S)]
        H = StS \ StV

        FShurr = schur(H)
        S = [sum(FShurr.vectors[:, i] .* V) for i in 1:dp.eigN]
        S .= S ./ norm.(S)
    end
    
    Eigvals = eigvals(H)
    p_idx = sortperm(Eigvals, by=abs, rev=true)
    return Eigvals[p_idx], S[p_idx]
end

# Utilities for ForwardDiff compatibility
function randsimilar(x::AbstractArray, N::Int)::Vector{typeof(x)}
    [rand!(copy(x)) for _ in 1:N]
end

function randsimilar(x::StaticArray, N::Int)::Vector{typeof(x)}
    rand(typeof(x), N)
end

function partialpart(xSA::SVector)
    bb = [x.partials[1] for x in xSA]
    return SVector(bb...)
end

function valuepart(xSA::SVector)
    bb = [x.value for x in xSA]
    return MVector(bb...)
end

function partialpart(xSA::StaticArray)
    bb = [x.partials[1] for x in xSA]
    return MVector(bb...)
end

function valuepart(xSA::StaticArray)
    bb = [x.value for x in xSA]
    return MVector(bb...)
end

function partialpart(xSA::AbstractArray)
    return [x.partials[1] for x in xSA]
end

function valuepart(xSA::AbstractArray)
    return [x.value for x in xSA]
end
