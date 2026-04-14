# AffineMapDiffEq.jl Mapping Logic

struct AffineTag end

"""
    dynamic_problemSampled(prob, alg, maxdelay; Historyresolution=200, zerofixpont=true, affineinteration=1, Krylov_arg=())

Constructor for `dynamic_problemSampled`. 
"""
function dynamic_problemSampled(prob, alg, maxdelay; 
                               Historyresolution=200,
                               zerofixpont=true, 
                               affineinteration=1, 
                               Krylov_arg=())
    
    if Historyresolution == 1
        StateSmaplingTime = [0.0]
    else
        StateSmaplingTime = LinRange(-maxdelay, 0.0, Historyresolution)
    end
    
    dynamic_problemSampled(prob, alg, maxdelay, zerofixpont, affineinteration,
        StateSmaplingTime, Krylov_arg)
end

"""
    LinMap(dp::dynamic_problemSampled, s::T; p=dp.Problem.p) where {T}

Core linear mapping function that integrates the DDE problem over one period.
"""
function LinMap(dp::dynamic_problemSampled, s::T; p=dp.Problem.p)::Tuple{T, ODESolution} where {T}
    StateSmaplingTime = dp.StateSmaplingTime
    
    # Use the type-stable history constructor
    ah = construct_history(s, StateSmaplingTime)

    # Standardized parameter passing and type-stable history wrapper
    new_prob = remake(dp.Problem; u0=ah(p, 0.0), h=ah, p=p)
    sol = solve(new_prob; dp.alg...)
    
    tend = sol.t[end]
    v = [getvalues(sol, ti + tend, p) for ti in StateSmaplingTime]
    
    return v, sol
end

"""
    affine(dp::dynamic_problemSampled, s0::T; p=dp.Problem.p, ...) where {T}

Performs affine mapping iteration to find the fixed point of the mapping.
"""
function affine(dp::dynamic_problemSampled, s0::T; 
                p=dp.Problem.p, 
                pDual_dir=nothing, 
                Δu=nothing, 
                Δλ_scaled=1.0, 
                norm_limit=1e-10,
                max_fix_iter=30) where {T}

    # Standardized ForwardDiff strategy with unique Tag
    one_epsilon_Dual = ForwardDiff.Dual{AffineTag}(0.0, 1.0)
    
    # Initialize directions if not provided
    p_dir = isnothing(pDual_dir) ? map(x -> 0.0, p) : pDual_dir
    du_dir = isnothing(Δu) ? map(x -> zero(x), s0) : Δu

    # Initial mapping
    v0, _ = LinMap(dp, s0; p=p)
    dv0dλ = map(x -> zero(x), v0) # Placeholder if not doing continuation

    Niteration = 0
    Finished_iteration = 0
    do_more_iteration = true
    
    while do_more_iteration
        Nstep = length(dp.StateSmaplingTime)
        s_start = randsimilar(dp.Problem.u0, Nstep)

        # Function barrier for the mapping perturbation
        TheMapping = let dp=dp, s0=s0, v0=v0, tag=AffineTag()
            function (ds::T) where {T}
                one_eps = ForwardDiff.Dual{typeof(tag)}(0.0, 1.0)
                perturbed_v, _ = LinMap(dp, ds .* one_eps .+ s0; p=p)
                return partialpart.(perturbed_v .- v0)
            end
        end

        # Spectral analysis
        mus = getindex(schursolve(TheMapping, s_start, dp.Krylov_arg...), [3, 2, 1])
        eigval, eigvec = mus[1], mus[2]

        # Fixed-point iteration
        a0, Δλ_loc = find_fix_pont(s0, v0, eigval, eigvec, dv0dλ, du_dir, Δλ_scaled)

        for _ in 1:max_fix_iter
            Niteration += 1
            s0 = a0
            v0, _ = LinMap(dp, s0; p=p)
            a0, Δλ_loc = find_fix_pont(s0, v0, eigval, eigvec, dv0dλ, du_dir, Δλ_scaled)
            if norm(s0 .- v0) < norm_limit
                break
            end
        end
        
        v0, sol = LinMap(dp, s0; p=p)
        a0, Δλ_loc = find_fix_pont(s0, v0, mus[1], mus[2], dv0dλ, du_dir, Δλ_scaled)

        s0 = a0
        v0, sol = LinMap(dp, s0; p=p)
        norm_err = norm(s0 .- v0)
        
        Finished_iteration += 1
        do_more_iteration = Finished_iteration < dp.affineinteration
        
        if !do_more_iteration
            return mus, s0::T, sol, p, Niteration, norm_err
        end
    end
end

function affine(dp::dynamic_problemSampled; p=dp.Problem.p, kwargs...)
    Nstep = length(dp.StateSmaplingTime)
    s0 = randsimilar(dp.Problem.u0, Nstep)
    if dp.zerofixpont
        fill!(s0, zero(eltype(s0)))
    end
    return affine(dp, s0; p=p, kwargs...)
end

# Helper functions
function find_fix_pont(s0::T, v0::T, eigval, eigvec) where {T}
    x = v0 .- s0
    AtA = [dot(eigvec[i], eigvec[j]) for i in eachindex(eigvec), j in eachindex(eigvec)]
    Atx = [dot(eigvec[i], x) for i in eachindex(eigvec)]
    ci = AtA \ Atx
    ci_mu = ci .* (eigval ./ (eigval .- 1.0))
    fix_v = v0 .- mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
    return fix_v, 0.0
end

function find_fix_pont(s0::T, v0::T, eigval, eigvec, dv0dλ, Δu, Δλ_scaled::Tlam) where {T,Tlam}
    x = v0 .- s0
    Atx = [dot(eigvec[i], x) for i in eachindex(eigvec)]
    Atdvdlam = [dot(eigvec[i], dv0dλ) for i in eachindex(eigvec)]
    ΔuAt = [dot(Δu, eigvec[i]) for i in eachindex(eigvec)]
    
    T_jac = vcat(hcat(diagm(eigval .- 1.0), Atdvdlam),
                 hcat(ΔuAt', Δλ_scaled))
    
    ci_arch = T_jac \ vcat(-Atx, 0.0)
    ci = ci_arch[1:end-1]
    Δλ_loc = real(ci_arch[end])
    ci_mu = ci .* eigval
    
    fix_v = real.(v0 .+ mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu)))
    return fix_v::T, Δλ_loc::Tlam
end

@inline function getvalues(sol::ODESolution, t, p)
    if t < 0.0
        return sol.prob.h(p, t)
    end
    return sol(t)
end
