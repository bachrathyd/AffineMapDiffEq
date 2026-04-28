# AffineMapDiffEq.jl Mapping Logic

struct AffineTag end

"""
    dynamic_problemSampled(prob, alg, maxdelay; Historyresolution=nothing, zerofixpont=true, affineinteration=1, Krylov_arg=(), perturbation_size=0.0, itp_type=DEFAULT_ITP)

Constructor for the `dynamic_problemSampled` configuration structure.

# Arguments
- `prob`: The underlying `DDEProblem` or `ODEProblem`.
- `alg`: Solver options as a `Dict`.
- `maxdelay`: The maximum delay (used for DDE history length). For ODEs, this can be 0.0.

# Keywords
- `Historyresolution`: Number of points used to discretize the history. Defaults to 1 for ODEs and 100 for DDEs.
- `zerofixpont`: If `true`, assumes the fixed point of the mapping is at the origin.
- `affineinteration`: Number of outer iterations for the affine mapping (re-calculating the spectrum).                                                                                                                                                    │    
- `Krylov_arg`: Arguments for the spectral solver (e.g., `(Neig, :LM, Arnoldi())`).                                                                                                                                                                       │    
- `perturbation_size`: Step size for finite difference mapping fallback. If `0.0`, `ForwardDiff` is used.    
- `itp_type`: Interpolation type from `Interpolations.jl`.
"""
function dynamic_problemSampled(prob, alg, maxdelay; 
                               Historyresolution=nothing,
                               zerofixpont=true, 
                               affineinteration=1, 
                               Krylov_arg=(),
                               perturbation_size=0.0,
                               itp_type=DEFAULT_ITP)
    
    is_ode = prob isa ODEProblem
    
    # Set default Historyresolution based on problem type
    if isnothing(Historyresolution)
        res = is_ode ? 1 : 100
    else
        res = Historyresolution
    end

    # Handle ODE specific constraints
    if is_ode && res > 1
        @warn "Problem is an ODEProblem but Historyresolution > 1. Setting Historyresolution = 1 for efficiency."
        res = 1
    end
    
    if res == 1
        StateSmaplingTime = [0.0]
    else
        StateSmaplingTime = LinRange(-maxdelay, 0.0, res)
    end
    
    dynamic_problemSampled(prob, alg, Float64(maxdelay), zerofixpont, affineinteration,
        StateSmaplingTime, Krylov_arg, Float64(perturbation_size), itp_type)
end

"""
    LinMap(dp::dynamic_problemSampled, s::T; p=dp.Problem.p) where {T}

Linear mapping operator \$\\mathcal{A}\$ that maps an initial history segment `s` to the history segment after one period \$T\$.

# Returns
- `v`: The resulting discretized history segment.
- `sol`: The `ODESolution` object for the integration period.
"""
function LinMap(dp::dynamic_problemSampled, s::T; p=dp.Problem.p)::Tuple{T, ODESolution} where {T}
    StateSmaplingTime = dp.StateSmaplingTime
    T_period = dp.Problem.tspan[2]
    
    # Check problem type
    is_ode = dp.Problem isa ODEProblem
    
    if is_ode
        # For ODEs, s should be a vector of length 1 containing u0
        new_prob = remake(dp.Problem; u0=s[1], p=p)
    else
        # Use the type-stable history constructor with configured interpolation for DDEs
        ah = construct_history(s, StateSmaplingTime; itp_type=dp.itp_type)
        # Standardized parameter passing and type-stable history wrapper
        new_prob = remake(dp.Problem; u0=ah(p, 0.0), h=ah, p=p)
    end
    
    # If the simulation period is longer than the delay, the required history 
    # for the next step is fully contained within the current simulation [0, T].
    # For ODEs, T_period >= 0.0 is always true since maxdelay is 0.0.
    if T_period >= dp.maxdelay
        save_times = StateSmaplingTime .+ T_period
        
        # Overwrite/add saveat and force save_everystep=false for maximum efficiency.
        sol = solve(new_prob; dp.alg..., saveat=save_times, save_everystep=false,dense=false)
        
        tend = sol.t[end]
        if tend >= T_period && length(sol.u) >= length(save_times)
            # Optimized extraction: 
           v = sol.u[1:end]
        else
            # Fallback for early exit or unexpected solver behavior
            v = map(ti -> getvalues(sol, ti + tend, p), StateSmaplingTime)
        end
    else
        # Period is shorter than delay: we need both the simulation result 
        # and the previous history (handled by getvalues).
        sol = solve(new_prob; dp.alg...)
        tend = sol.t[end]
        v = map(ti -> getvalues(sol, ti + tend, p), StateSmaplingTime)
    end
    
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
                max_fix_iter=30,
                method=:Krylov) where {T}

    # Initialize directions if not provided
    du_dir = isnothing(Δu) ? map(x -> zero(x), s0) : Δu

    # Initial mapping
    v0, _ = LinMap(dp, s0; p=p)
    dv0dλ = map(x -> zero(x), v0) # Placeholder if not doing continuation

    Niteration = 0
    Finished_iteration = 0
    do_more_iteration = true
    
    mus = nothing
    sol = v0 # Placeholder for sol if loop doesn't run
    norm_err = norm(s0 .- v0)

    while do_more_iteration
        Nstep = length(dp.StateSmaplingTime)
        s_start = randsimilar(dp.Problem.u0, Nstep)

        # Function barrier for the mapping perturbation
        TheMapping = let dp=dp, s0=s0, tag=AffineTag(), p=p
            function (ds::T) where {T}
                if dp.perturbation_size == 0.0
                    # Exact Jacobian action via ForwardDiff directional derivative
                    one_eps = ForwardDiff.Dual{typeof(tag)}(0.0, 1.0)
                    # Use efficient broadcasting for perturbation
                    perturbed_s = s0 .+ ds .* one_eps
                    perturbed_v, _ = LinMap(dp, perturbed_s; p=p)
                    # Extract directional derivative
                    return map(partialpart, perturbed_v)
                else
                    # Finite difference fallback
                    eps_val = dp.perturbation_size
                    perturbed_v, _ = LinMap(dp, s0 .+ ds .* eps_val; p=p)
                    v_base, _ = LinMap(dp, s0; p=p)
                    return (perturbed_v .- v_base) ./ eps_val
                end
            end
        end

        # Spectral analysis
        if method == :Krylov
            eigval, eigvec, info = spectrum_krylov(TheMapping, s_start, dp.Krylov_arg...)
        elseif method == :ISSI
            Neig = dp.Krylov_arg[1]
            eigval, eigvec, info = issi_eigen(TheMapping, s_start, Neig)
        else
            error("Unknown method: $method")
        end
        mus = (eigval, eigvec, info)

        # Fixed-point iteration
        a0, Δλ_loc = find_fixed_point(s0, v0, eigval, eigvec, dv0dλ, du_dir, Δλ_scaled)

        for _ in 1:max_fix_iter
            Niteration += 1
            s0 = a0
            v0, sol = LinMap(dp, s0; p=p) # Update sol here
            a0, Δλ_loc = find_fixed_point(s0, v0, eigval, eigvec, dv0dλ, du_dir, Δλ_scaled)
            norm_err = norm(s0 .- v0)
            if norm_err < norm_limit
                break
            end
        end
        
        
        # Removed redundant LinMap calls that were here
        
        Finished_iteration += 1
        do_more_iteration = Finished_iteration < dp.affineinteration
    end
    
    return mus, s0, sol, p, Niteration, norm_err
end

function affine(dp::dynamic_problemSampled; p=dp.Problem.p, kwargs...)
    Nstep = length(dp.StateSmaplingTime)
    s0 = [zero(dp.Problem.u0) for _ in 1:Nstep]
    if !dp.zerofixpont
        # Initialize with random if not zero fixed point and not provided
        s0 = randsimilar(dp.Problem.u0, Nstep)
    end
    return affine(dp, s0; p=p, kwargs...)
end

# Helper functions
function find_fixed_point(s0::T, v0::T, eigval, eigvec) where {T}
    x = v0 .- s0
    # Use real part if eigenvalues are complex but result should be real? 
    # Actually keep it complex if needed, but usually we want real fixed point.
    AtA = [dot(eigvec[i], eigvec[j]) for i in eachindex(eigvec), j in eachindex(eigvec)]
    Atx = [dot(eigvec[i], x) for i in eachindex(eigvec)]
    ci = AtA \ Atx
    ci_mu = ci .* (eigval ./ (eigval .- 1.0))
    
    # Efficiently sum the corrections
    correction = mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
    fix_v = v0 .- correction
    
    return fix_v, 0.0
end

function find_fixed_point(s0::T, v0::T, eigval, eigvec, dv0dλ, Δu, Δλ_scaled::Tlam) where {T,Tlam}
    x = v0 .- s0
    Neig = length(eigval)
    
    # Pre-calculate dot products to avoid multiple passes
    Atx = [dot(eigvec[i], x) for i in 1:Neig]
    Atdvdlam = [dot(eigvec[i], dv0dλ) for i in 1:Neig]
    ΔuAt = [dot(Δu, eigvec[i]) for i in 1:Neig]
    
    # Construct Jacobian matrix for the fixed point / continuation problem
    # [ diag(eigval - 1)  A*dvdlam ]
    # [ ΔuAt'              Δλ_scaled ]
    T_jac = zeros(ComplexF64, Neig + 1, Neig + 1)
    for i in 1:Neig
        T_jac[i, i] = eigval[i] - 1.0
        T_jac[i, Neig+1] = Atdvdlam[i]
        T_jac[Neig+1, i] = ΔuAt[i]
    end
    T_jac[Neig+1, Neig+1] = Δλ_scaled
    
    # RHS: [-Atx, 0]
    rhs = zeros(ComplexF64, Neig + 1)
    for i in 1:Neig
        rhs[i] = -Atx[i]
    end
    
    ci_arch = T_jac \ rhs
    ci = ci_arch[1:Neig]
    Δλ_loc = real(ci_arch[end])
    ci_mu = ci .* eigval
    
    # Summing up the corrections
    correction = mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
    fix_v = v0 .+ correction
    
    # Return real part if the state type is Real
    if eltype(eltype(s0)) <: Real
        return real.(fix_v), Δλ_loc
    else
        return fix_v, Δλ_loc
    end
end

@inline function getvalues(sol::ODESolution, t, p)
    if t < 0.0
        # For DDEs, sol.prob.h is the history function
        # Check if it's a DDEProblem by checking if .h exists
        if hasfield(typeof(sol.prob), :h)
            return sol.prob.h(p, t)
        else
            # For ODEs, if t < 0 it's technically undefined but we can return u0 
            # or handle it based on how StateSmaplingTime is constructed.
            return sol.prob.u0
        end
    end
    return sol(t)
end
