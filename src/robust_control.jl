# Robustness margin C(p) and controller optimisation for Floquet-based stability
#
#   C_j = (1 - |μ_j(p₀)|) / Σ_i |∂|μ_j|/∂p_i| Δp_i     (per multiplier)
#   C   = min_j C_j                                         (worst-case)
#
# Sensitivity ∂|μ_j|/∂p_i via central FD + root tracking.
# The nominal fixed-point history is reused as a warm start for every perturbed
# affine call — critical for DDE systems where zerofixpont=false.

using LinearAlgebra
using ForwardDiff

# ─── Root Tracking ────────────────────────────────────────────────────────────

"""
    track_roots(mu_ref, mu_new) -> Vector

Greedy O(n²) bipartite matching of `mu_new` onto `mu_ref` by minimum |Δμ|.
Prevents spurious eigenvalue branch crossings under small parameter perturbations.
"""
function track_roots(mu_ref::AbstractVector{<:Complex},
                     mu_new::AbstractVector{S})::Vector{S} where {S<:Complex}
    n = length(mu_ref)
    @assert length(mu_new) == n "Root count mismatch: $(length(mu_new)) ≠ $n"
    result = Vector{S}(undef, n)
    used   = falses(n)
    @inbounds for i in 1:n
        best_j, best_d2 = 0, Inf
        for j in 1:n
            used[j] && continue
            d2 = abs2(mu_ref[i] - mu_new[j])
            if d2 < best_d2
                best_d2 = d2
                best_j  = j
            end
        end
        result[i]    = mu_new[best_j]
        used[best_j] = true
    end
    return result
end

# ─── Eigenvalue Extraction ────────────────────────────────────────────────────

"""
    get_eigenvalues(dp, p, n_mu, s_warmstart) -> Vector{ComplexF64}

Run `affine` at parameter tuple `p`, starting from `s_warmstart` as initial
history estimate (warm start). Returns the `n_mu` dominant Floquet multipliers.

Using the nominal fixed point as warm start is essential for DDE systems
(zerofixpont=false): it avoids divergence from random initialisations and
reduces the fixed-point iteration count to 1–2 steps.
"""
function get_eigenvalues(dp, p::Tuple, n_mu::Int, s_warmstart)::Vector{ComplexF64}
    mus, = affine(dp, s_warmstart; p=p)   # mus = (eigvals, eigvecs, info)
    ev = ComplexF64.(mus[1])
    # Krylov may converge fewer than n_mu multipliers (common in the noisy FD
    # branch); return whatever is available rather than indexing out of bounds.
    return ev[1:min(n_mu, length(ev))]
end

# ─── Parameter Sensitivity ────────────────────────────────────────────────────
#
# Two branches, selected by `dp.perturbation_size` (mirrors the Jacobian-action
# logic in mapping.jl):
#   - perturbation_size != 0  →  central finite differences  (compute_sensitivities_fd)
#   - perturbation_size == 0  →  ForwardDiff Dual numbers     (compute_sensitivities_dual)
#
# The Dual branch differentiates the period map's compressed operator w.r.t. the
# parameters, so the sensitivities carry no FD step error or cancellation noise
# from the slightly inaccurate Krylov |μ|. Accuracy improves with the Schur
# subspace size (Neig); see compute_sensitivities_dual.

# Tag for the outer (parameter) Dual layer. Distinct from `AffineTag` used inside
# the mapping for the Jacobian action, so the two Dual layers never mix.
struct SensTag end

# Nesting order for the two Dual layers. `partialpart` (inside the mapping) strips
# the OUTERMOST tag to extract the Jacobian-action direction, so `AffineTag` must
# be the outer layer and the parameter-sensitivity `SensTag` the inner one — then
# stripping AffineTag leaves the SensTag (parameter) derivative intact. ForwardDiff
# needs this precedence to combine the two tags without an "ordering" error.
ForwardDiff.:≺(::Type{SensTag}, ::Type{AffineMapDiffEq.AffineTag}) = true
ForwardDiff.:≺(::Type{AffineMapDiffEq.AffineTag}, ::Type{SensTag}) = false

# Action of the linearized period map M = ∂F/∂s at the (Float64) fixed point
# `s_lin` on a (Float64) direction `ds`, with Dual-valued parameters `p_dual`.
# Mirrors the inner-AffineTag Dual mechanism of `TheMapping` in mapping.jl: the
# AffineTag directional derivative gives M·ds, and the SensTag partials it carries
# give ∂(M·ds)/∂p_i. Returns a state vector whose elements are Dual{SensTag}.
#
# `one_eps` is the AffineTag unit-direction seed, but nested as
# `Dual{AffineTag, Dual{SensTag,Float64,N}, 1}` with zero SensTag partials. This
# makes the constructed history `Dual{SensTag}`-wide, so when a parameter that
# enters the delayed time argument (e.g. τ) is itself a `Dual{SensTag}`, querying
# the history at the resulting Dual time stays within the declared element type.
function _linmap_action_dual(dp, s_lin, ds, p_dual, one_eps)
    perturbed_s    = s_lin .+ ds .* one_eps
    perturbed_v, _ = LinMap(dp, perturbed_s; p=p_dual)
    return map(partialpart, perturbed_v)
end

"""
    compute_sensitivities(dp, p_nominal, delta_p, mu_nominal, s_nominal, eigvecs; delta_rel)
        -> Matrix{Float64}

Return `sens[j, i] = ∂|μ_j|/∂p_i`. Dispatches on `dp.perturbation_size`:
exact ForwardDiff (`== 0.0`) or central finite differences (`!= 0.0`).
`eigvecs` is the nominal Schur basis (only used by the ForwardDiff branch).
Columns where `delta_p[i] == 0` are skipped.
`verbose=false` suppresses progress prints (required inside a threaded sweep).
"""
function compute_sensitivities(dp,
                                p_nominal::Tuple,
                                delta_p::Tuple,
                                mu_nominal::AbstractVector{<:Complex},
                                s_nominal,
                                eigvecs;
                                delta_rel::Float64 = 1e-5,
                                verbose::Bool = true)::Matrix{Float64}
    if dp.perturbation_size == 0.0
        return compute_sensitivities_dual(dp, p_nominal, delta_p, mu_nominal, s_nominal, eigvecs)
    else
        return compute_sensitivities_fd(dp, p_nominal, delta_p, mu_nominal, s_nominal; delta_rel, verbose)
    end
end

"""
    compute_sensitivities_fd(dp, p_nominal, delta_p, mu_nominal, s_nominal; delta_rel, verbose)

2nd-order central finite differences + root tracking.
FD step `δ = max(|p_i|, 1) × delta_rel`. Cost: 2 × (active parameters) warm-started `affine` calls.
"""
function compute_sensitivities_fd(dp,
                                   p_nominal::Tuple,
                                   delta_p::Tuple,
                                   mu_nominal::AbstractVector{<:Complex},
                                   s_nominal;
                                   delta_rel::Float64 = 1e-5,
                                   verbose::Bool = true)::Matrix{Float64}
    n_params = length(p_nominal)
    n_mu     = length(mu_nominal)
    p_vec    = collect(Float64, p_nominal)
    sens     = zeros(Float64, n_mu, n_params)

    for i in 1:n_params
        iszero(delta_p[i]) && continue

        δ     = max(abs(p_vec[i]), 1.0) * delta_rel
        p_fwd = ntuple(k -> k == i ? p_vec[k] + δ : p_vec[k], n_params)
        p_bwd = ntuple(k -> k == i ? p_vec[k] - δ : p_vec[k], n_params)

        ev_fwd = get_eigenvalues(dp, p_fwd, n_mu, s_nominal)
        ev_bwd = get_eigenvalues(dp, p_bwd, n_mu, s_nominal)

        # The perturbed Krylov solves may return fewer multipliers than n_mu;
        # only difference the modes available in both, leaving the rest at 0.
        nf = min(length(ev_fwd), length(ev_bwd), n_mu)
        nf == 0 && continue
        mu_fwd = track_roots(mu_nominal[1:nf], ev_fwd[1:nf])
        mu_bwd = track_roots(mu_nominal[1:nf], ev_bwd[1:nf])

        @inbounds for j in 1:nf
            sens[j, i] = (abs(mu_fwd[j]) - abs(mu_bwd[j])) / (2δ)
        end

        verbose && print("  param $i/$n_params done\r")
    end
    verbose && println()
    return sens
end

"""
    compute_sensitivities_dual(dp, p_nominal, delta_p, mu_nominal, s_nominal, eigvecs)

`∂|μ_j|/∂p_i` via ForwardDiff, without ever running the Krylov solver on Dual
numbers (`KrylovKit.schursolve` only supports BLAS floats).

Strategy — valid because the systems here are linear in the state, so the
monodromy operator `M` depends on `p` only through the coefficients, not through
the fixed-point location:

1. Take the orthonormal Schur basis `Q = eigvecs` from the nominal Float64 solve.
2. Seed `p` as chunked `Dual{SensTag}` numbers (one partial per active parameter).
3. Build the compressed operator `H[a,b] = ⟨Q_a, M(p_dual)·Q_b⟩` with `nb` Dual
   `LinMap` calls. `H` is a tiny `nb×nb` matrix carrying parameter partials.
4. Read `∂|λ_m|/∂p_i` off `H` by first-order eigenvalue perturbation of the
   small matrix — all in Float64 linear algebra (`eigen` on `value.(H)`).

The result is free of FD step/cancellation noise. It is first-order exact in the
limit where `Q` spans the full left/right eigenspaces; for a truncated subspace a
small bias remains in `∂λ` because the modes' left eigenvectors are only captured
inside `span(Q)`. Increasing the subspace size `Neig` shrinks this bias quickly
(empirically near-exact by `Neig ≈ 10–12` for the 2-DOF examples here).

Cost: `nb` warm Dual `LinMap` calls (vs 2×active full `affine` calls for FD).
"""
function compute_sensitivities_dual(dp,
                                     p_nominal::Tuple,
                                     delta_p::Tuple,
                                     mu_nominal::AbstractVector{<:Complex},
                                     s_nominal,
                                     eigvecs)::Matrix{Float64}
    n_params = length(p_nominal)
    n_mu     = length(mu_nominal)
    p_vec    = collect(Float64, p_nominal)
    sens     = zeros(Float64, n_mu, n_params)

    active   = Int[i for i in 1:n_params if !iszero(delta_p[i])]
    n_active = length(active)
    n_active == 0 && return sens

    DualT = ForwardDiff.Dual{SensTag, Float64, n_active}
    zp    = ForwardDiff.Partials(ntuple(_ -> 0.0, n_active))
    p_dual = ntuple(n_params) do k
        idx = findfirst(==(k), active)
        if idx === nothing
            DualT(p_vec[k], zp)
        else
            DualT(p_vec[k], ForwardDiff.Partials(ntuple(j -> j == idx ? 1.0 : 0.0, n_active)))
        end
    end

    # AffineTag unit-direction seed, nested over Dual{SensTag} with zero partials.
    one_eps = ForwardDiff.Dual{AffineMapDiffEq.AffineTag, DualT, 1}(
        DualT(0.0, zp), ForwardDiff.Partials((DualT(1.0, zp),)))

    # Compressed operator H[a,b] = ⟨Q_a, M(p_dual) Q_b⟩ (Dual-valued, nb×nb)
    Q  = eigvecs
    nb = length(Q)
    MQ = [_linmap_action_dual(dp, s_nominal, Q[b], p_dual, one_eps) for b in 1:nb]
    Hd = [dot(Q[a], MQ[b]) for a in 1:nb, b in 1:nb]

    # Value part → small eigenproblem; left/right eigenvectors for perturbation
    H0  = ForwardDiff.value.(Hd)
    F   = eigen(H0)
    λ   = F.values
    R   = F.vectors
    Lt  = inv(R)                       # rows = left eigenvectors, Lt*R = I

    # ∂|λ_m|/∂p_i for each computed mode m and active parameter i
    dabs = zeros(Float64, nb, n_params)
    for (k_idx, i) in enumerate(active)
        H1 = ForwardDiff.partials.(Hd, k_idx)      # ∂H/∂p_i  (nb×nb, Float64)
        M1 = Lt * H1 * R                            # ∂λ in the eigenbasis (diagonal = ∂λ_m)
        for m in 1:nb
            am = abs(λ[m])
            dabs[m, i] = am == 0.0 ? 0.0 : real(conj(λ[m]) * M1[m, m]) / am
        end
    end

    # Match computed modes λ onto the nominal multiplier ordering
    used = falses(nb)
    for j in 1:n_mu
        best_m, best_d2 = 0, Inf
        for m in 1:nb
            used[m] && continue
            d2 = abs2(mu_nominal[j] - λ[m])
            if d2 < best_d2
                best_d2, best_m = d2, m
            end
        end
        best_m == 0 && continue
        used[best_m] = true
        @inbounds for i in active
            sens[j, i] = dabs[best_m, i]
        end
    end

    return sens
end

# ─── Robustness Margin ────────────────────────────────────────────────────────

"""
    robustness_margin(dp, p_nominal, delta_p; n_mu, delta_rel)
        -> (C, C_vec, mu_nominal, sens)

Worst-case linear robustness margin under bounded parameter uncertainty.

  C_j = (1 - |μ_j|) / Σ_i |∂|μ_j|/∂p_i| · Δp_i
  C   = min_j C_j

Returns `(-Inf, ...)` if any |μ_j| ≥ 1 (unstable nominal point).
Returns `Inf` for a root whose sensitivity denominator is numerically zero.

The nominal `affine` call produces both the Floquet multipliers and the
periodic-orbit fixed point, which is passed as a warm start to all
subsequent perturbed evaluations.
"""
function robustness_margin(dp,
                            p_nominal::Tuple,
                            delta_p::Tuple;
                            n_mu::Int          = dp.Krylov_arg[1],
                            delta_rel::Float64 = 1e-5,
                            verbose::Bool      = true)
    # ── Nominal computation (also produces warm-start history s_nominal) ──
    mus_nominal, s_nominal = affine(dp; p=p_nominal)
    all_mu     = ComplexF64.(mus_nominal[1])
    n_eff      = min(n_mu, length(all_mu))   # Krylov may converge fewer than n_mu
    mu_nominal = all_mu[1:n_eff]
    eigvecs    = mus_nominal[2]
    mu_abs     = abs.(mu_nominal)

    if any(>=(1.0), mu_abs)
        sens_empty = zeros(Float64, n_eff, length(p_nominal))
        return -Inf, fill(-Inf, n_eff), mu_nominal, sens_empty
    end

    # ── Sensitivities using s_nominal as warm start ──
    if verbose
        if dp.perturbation_size == 0.0
            println("  Computing sensitivities (ForwardDiff on compressed operator)...")
        else
            println("  Computing sensitivities (central FD, $(2 * count(!iszero, delta_p)) affine calls)...")
        end
    end
    sens      = compute_sensitivities(dp, p_nominal, delta_p, mu_nominal, s_nominal, eigvecs; delta_rel, verbose)
    delta_vec = collect(Float64, delta_p)

    C_vec = map(eachindex(mu_nominal)) do j
        denom = sum(abs(sens[j, i]) * delta_vec[i] for i in eachindex(delta_vec))
        denom < eps(Float64) ? Inf : (1.0 - mu_abs[j]) / denom
    end

    return minimum(C_vec), C_vec, mu_nominal, sens
end

# TODO: finalize it for controller optimization (Optim pkg is needed)
#
# """
#     optimize_controller(dp, p_base, delta_p, ctrl_indices;
#                         k_bounds, n_mu, method, optim_options)
#         -> (k_opt, C_opt, opt_result)
#
# Maximize C(k) over controller gains at `ctrl_indices` inside the parameter tuple.
#
# - `ctrl_indices`: e.g. `[1, 4]` for gains at positions 1 and 4 of the tuple
# - `k_bounds`    : `(k_lo, k_hi)` vectors — required for `:ParticleSwarm`
# - `method`      : `:NelderMead` (default) | `:ParticleSwarm`
# """
# function optimize_controller(dp,
#                               p_base::Tuple,
#                               delta_p::Tuple,
#                               ctrl_indices::AbstractVector{Int};
#                               k_bounds      = nothing,
#                               n_mu::Int     = dp.Krylov_arg[1],
#                               method        = :NelderMead,
#                               optim_options = Optim.Options(iterations = 300,
#                                                             f_tol      = 1e-4,
#                                                             show_trace = true))
#     p_vec  = collect(Float64, p_base)
#     n_ctrl = length(ctrl_indices)
#     k0     = p_vec[ctrl_indices]
#
#     function neg_C(k_vec::Vector{Float64})::Float64
#         p_eval = copy(p_vec)
#         @inbounds p_eval[ctrl_indices] .= k_vec
#         C, = robustness_margin(dp, Tuple(p_eval), delta_p; n_mu)
#         return -C
#     end
#
#     if method == :NelderMead
#         if isnothing(k_bounds)
#             opt_result = optimize(neg_C, k0, NelderMead(), optim_options)
#         else
#             k_lo, k_hi = k_bounds
#             function penalized(k_vec::Vector{Float64})::Float64
#                 pen = 1e6 * sum(max(0.0, k_lo[i] - k_vec[i])^2 +
#                                 max(0.0, k_vec[i] - k_hi[i])^2 for i in 1:n_ctrl)
#                 return neg_C(k_vec) + pen
#             end
#             opt_result = optimize(penalized, k0, NelderMead(), optim_options)
#         end
#
#     elseif method == :ParticleSwarm
#         isnothing(k_bounds) && error("ParticleSwarm requires k_bounds = (k_lo, k_hi)")
#         k_lo, k_hi = k_bounds
#         ps_alg     = ParticleSwarm(lower = k_lo, upper = k_hi, n_particles = 20)
#         opt_result = optimize(neg_C, k_lo, k_hi, k0, ps_alg, optim_options)
#
#     else
#         error("Unknown method :$method — choose :NelderMead or :ParticleSwarm")
#     end
#
#     k_opt = Optim.minimizer(opt_result)
#     C_opt = -Optim.minimum(opt_result)
#     return k_opt, C_opt, opt_result
# end
