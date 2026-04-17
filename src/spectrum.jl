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
Optimized to avoid unnecessary allocations.
"""
function issi_eigen(f, s_start::T, Neig; max_iter=12) where {T}
    Nstep = length(s_start)
    s0 = map(x -> zero(x), s_start)
    v0 = f(s0)
    
    S = [randsimilar(s0[1], Nstep) for _ in 1:Neig]
    V = [similar(S[1]) for _ in 1:Neig]
    
    for _ in 1:max_iter
        for i in 1:Neig
            V[i] .= f(S[i] .+ s0) .- v0
        end
        
        # Compute H = (S'S) \ (S'V) efficiently using dot products
        Gram = [dot(S[i], S[j]) for i in 1:Neig, j in 1:Neig]
        SV = [dot(S[i], V[j]) for i in 1:Neig, j in 1:Neig]
        
        H = Gram \ SV

        FShurr = schur(H)
        for i in 1:Neig
            # Recompute S as a linear combination of V using Schur vectors
            S[i] .= sum(FShurr.vectors[j, i] .* V[j] for j in 1:Neig)
            S[i] ./= norm(S[i])
        end
    end
    
    # Final eigenvalue extraction from the projected subspace
    H_final = [dot(S[i], f(S[j] .+ s0) .- v0) for i in 1:Neig, j in 1:Neig]
    Eigvals, Eigvecs_proj = eigen(H_final)
    
    # Reconstruct eigenvectors in the original space
    S_final = [sum(Eigvecs_proj[j, i] .* S[j] for j in 1:Neig) for i in 1:Neig]
    
    p_idx = sortperm(Eigvals, by=abs, rev=true)
    return Eigvals[p_idx], S_final[p_idx], nothing
end

# Utilities for ForwardDiff compatibility
function randsimilar(x::AbstractArray, N::Int)::Vector{typeof(x)}
    [rand!(copy(x)) for _ in 1:N]
end

function randsimilar(x::StaticArray, N::Int)::Vector{typeof(x)}
    [rand(typeof(x)) for _ in 1:N]
end

function partialpart(xSA::StaticArray)
    return map(x -> ForwardDiff.partials(x, 1), xSA)
end

function valuepart(xSA::StaticArray)
    return map(x -> ForwardDiff.value(x), xSA)
end

function partialpart(xSA::AbstractArray)
    return [ForwardDiff.partials(x, 1) for x in xSA]
end

function valuepart(xSA::AbstractArray)
    return [ForwardDiff.value(x) for x in xSA]
end

function partialpart(x::ForwardDiff.Dual)
    return ForwardDiff.partials(x, 1)
end

function valuepart(x::ForwardDiff.Dual)
    return ForwardDiff.value(x)
end

function partialpart(x::Number)
    return zero(x)
end

function valuepart(x::Number)
    return x
end
