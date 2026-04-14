# AffineMapDiffEq.jl Spectrum and Eigen Analysis

"""
    spectrum(dp::dynamic_problemSampled; p=dp.Problem.p, method=:Krylov)

Calculates the Floquet multipliers (spectrum) of the mapping.

# Methods
- `:Krylov`: Uses `schursolve` from `KrylovKit.jl` (default).
- `:ISSI`: Uses Internal Subspace Iteration (ISSI).

# Returns
- `vals`: Eigenvalues (Floquet multipliers).
- `vecs`: Eigenvectors (discretized eigenfunctions).
- `info`: Convergence information from the solver.
"""
function spectrum(dp::dynamic_problemSampled; p=dp.Problem.p, method=:Krylov)
    Nstep = length(dp.StateSmaplingTime)
    s_start = randsimilar(dp.Problem.u0, Nstep)
    
    if method == :Krylov
        return spectrum_krylov(s -> LinMap(dp, s; p=p)[1], s_start, dp.Krylov_arg...)
    elseif method == :ISSI
        Neig = dp.Krylov_arg[1]
        return issi_eigen(s -> LinMap(dp, s; p=p)[1], s_start, Neig)
    else
        error("Unknown method: $method")
    end
end

"""
    spectrum_krylov(f, x0, args...)

Helper to call `schursolve` and extract eigenvalues and eigenvectors.
Returns `(vals, vecs, info)`.
"""
function spectrum_krylov(f, x0, args...)
    res = schursolve(f, x0, args...)
    # res = (T, vecs, vals, info)
    return res[3], res[2], res[4]
end

"""
    spectralradius(dp::dynamic_problemSampled; p=dp.Problem.p, method=:Krylov)

Calculates the maximum absolute Floquet multiplier, \$\\rho(\\mathcal{A})\$, which determines the stability of the periodic solution.
Stability is guaranteed if \$\\rho(\\mathcal{A}) < 1\$.
"""
function spectralradius(dp::dynamic_problemSampled; p=dp.Problem.p, method=:Krylov)
    if dp.zerofixpont
        return Float64(maximum(abs.(spectrum(dp; p=p, method=method)[1])))::Float64
    else
        return Float64(maximum(abs.(affine(dp; p=p)[1][1])))::Float64
    end
end

"""
    issi_eigen(f, s_start, Neig; max_iter=12)

Internal Subspace Iteration (ISSI) based eigen analysis.
"""
function issi_eigen(f, s_start::T, Neig; max_iter=12) where {T}
    Nstep = length(s_start)
    s0 = map(x -> zero(x), s_start)
    v0 = f(s0)
    
    S = [randsimilar(s0[1], Nstep) for _ in 1:Neig]
    H = zeros(ComplexF64, Neig, Neig)
    
    V = copy(S)
    for _ in 1:max_iter
        for i in 1:Neig
            V[i] = f(S[i] .+ s0) .- v0
        end
        
        # StS = [S[i]' * S[j] for i in 1:Neig, j in 1:Neig]
        # StV = [S[i]' * V[j] for i in 1:Neig, j in 1:Neig]
        # H = StS \ StV
        
        # More efficient way to compute H
        S_mat = reduce(hcat, S)
        V_mat = reduce(hcat, V)
        H = (S_mat' * S_mat) \ (S_mat' * V_mat)

        FShurr = schur(H)
        for i in 1:Neig
            S[i] = sum(FShurr.vectors[j, i] .* V[j] for j in 1:Neig)
            S[i] ./= norm(S[i])
        end
    end
    
    Eigvals = eigvals(H)
    p_idx = sortperm(Eigvals, by=abs, rev=true)
    return Eigvals[p_idx], S[p_idx], nothing # Return nothing for info
end

# Utilities for ForwardDiff compatibility
function randsimilar(x::AbstractArray, N::Int)::Vector{typeof(x)}
    [rand!(copy(x)) for _ in 1:N]
end

function randsimilar(x::StaticArray, N::Int)::Vector{typeof(x)}
    [rand(typeof(x)) for _ in 1:N]
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

function partialpart(x::ForwardDiff.Dual)
    return x.partials[1]
end

function valuepart(x::ForwardDiff.Dual)
    return x.value
end

function partialpart(x::Number)
    return zero(x)
end

function valuepart(x::Number)
    return x
end
