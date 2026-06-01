# Single-Point Sensitivity Analysis Demo
# Analyses the robustness of the delayed PD-controlled oscillator at ONE chosen
# controller setting (kp, kd), comparing the first-order (C₁) and second-order
# (C₂) robustness margins.
#
#   ẍ(t) + 2ζω ẋ(t) + ω² x(t) = -kp·x(t-τ) - kd·ẋ(t-τ) + F(t)
#
# Robustness margin C scales the uncertainty box so the worst-case |μ| reaches 1:
#   C₁ = (1 - |μ|) / Σ|∂|μ|/∂p|·Δp                (linear / gradient only)
#   C₂ : positive root of  ρ + C·L₁ + ½C²·L₂ = 1  (quadratic, Hessian, vertex check)
#
# Output: per-multiplier gradient, Hessian, and the C₁ / C₂ margins, plus a small
# Monte-Carlo check that perturbing the parameters by C·Δp brings |μ| to ≈1.
#
# Theory: docs/second_order_sensitivity_theory.md

using Pkg; Pkg.activate(".")
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using Statistics
using Printf

## ─── Model ────────────────────────────────────────────────────────────────────
function DelayPDOsc!(du, u, h, p, t)
    ζ, ω, kp, kd, τ, T = p
    hτ = h(p, t - τ)
    F  = 0.1 * cos(2π * t / T)
    du[1] = u[2]
    #Unstable with +ω^2, stable with -ω^2
    du[2] = +ω^2 * u[1] - 2ζ*ω * u[2] - kp * hτ[1] - kd * hτ[2] + F
end

const ζ_nom, ω_nom, τ_nom, T_nom = -0.02, 1.0, 0.5, 0.5

# ── The single controller point under analysis ────────────────────────────────
# Edit these two gains to analyse a different point of the (kp, kd) plane.
const kp_point = 1.2
const kd_point = 1.2

const p_point = (ζ_nom, ω_nom, kp_point, kd_point, τ_nom, T_nom)

u0 = @MArray [0.0, 0.0]
h(p, t) = @MArray [0.0, 0.0]
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-7)
probMapping = DDEProblem{true}(DelayPDOsc!, u0, h, (0.0, T_nom), p_point;
    constant_lags=[τ_nom])

include(joinpath(@__DIR__, "..", "src", "robust_control.jl"))

## ─── Uncertainty box (ratio of |nominal|) ──────────────────────────────────────
const ζ_var, ω_var, kp_var, kd_var, τ_var, T_var = 0.02, 0.00, 0.05, 0.05, 0.02, 0.0
const delta_p = (ζ_var*abs(ζ_nom), ω_var*abs(ω_nom), kp_var*abs(kp_point),
                 kd_var*abs(kd_point), τ_var*abs(τ_nom), 0.0)
const param_names = ("ζ", "ω", "kp", "kd", "τ", "T")

## ─── Dynamic problem (Dual sensitivities → needs perturbation_size = 0) ─────────
const Neig = 5
KA = (Neig, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig+20, verbosity=0, maxiter=20))
dp = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1,
    Krylov_arg=KA, perturbation_size=0.0)

## ─── Analysis at the single point ───────────────────────────────────────────────
println("="^78)
println(" Single-point sensitivity:  kp = $kp_point ,  kd = $kd_point")
println("="^78)
println("Uncertainty half-widths Δp:")
for i in eachindex(delta_p)
    @printf("   Δ%-2s = %.5f\n", param_names[i], delta_p[i])
end

# Warmup (JIT) then compute.
robustness_margin(dp, p_point, delta_p; verbose=false, sens_order=1)
robustness_margin(dp, p_point, delta_p; verbose=false, sens_order=2)

println("\nComputing first-order margin C₁...")
C1, C1_vec, mu, sens1 = robustness_margin(dp, p_point, delta_p; verbose=false, sens_order=1)

println("Computing second-order margin C₂...")
C2, C2_vec, _, _ = robustness_margin(dp, p_point, delta_p; verbose=false, sens_order=2)

# Raw gradient + Hessian of |μ| for inspection (per multiplier).
# affine returns (mus, s_nominal) with mus = (eigvals, eigvecs, info).
mus_pt, s_nom = affine(dp; p=p_point)
g_all, H_all, _ = compute_second_order_sens(dp, p_point, delta_p, mu, s_nom, mus_pt[2])

ρ = maximum(abs.(mu))
unstable = ρ >= 1.0

## ─── Report ──────────────────────────────────────────────────────────────────
println("\n── Floquet multipliers ────────────────────────────────────────────────────")
for j in eachindex(mu)
    @printf("   μ%-2d = %+.5f %+.5fi   |μ| = %.6f\n", j, real(mu[j]), imag(mu[j]), abs(mu[j]))
end
@printf("\n   ρ = max|μ| = %.6f   →   %s\n", ρ, unstable ? "UNSTABLE" : "stable")

if unstable
    println("\n   Nominal point is unstable; robustness margin is not defined (C = -Inf).")
else
    # Index of the worst (binding) multiplier for each order.
    j1 = argmin(C1_vec)
    j2 = argmin(C2_vec)

    println("\n── Sensitivity of the worst (binding) multiplier ──────────────────────────")
    @printf("   1st order binds on μ%-2d ;  2nd order binds on μ%-2d\n", j1, j2)
    println("\n   Gradient  ∂|μ$(j1)|/∂p_i :")
    for i in eachindex(delta_p)
        delta_p[i] == 0 && continue
        @printf("      ∂|μ|/∂%-2s = %+.5f\n", param_names[i], g_all[j1, i])
    end

    active = findall(!iszero, collect(delta_p))
    println("\n   Hessian  ∂²|μ$(j1)|/∂p_i∂p_k  (active params $(Tuple(param_names[i] for i in active))):")
    for i in active
        row = join((@sprintf("%+.4f", H_all[j1, i, k]) for k in active), "  ")
        @printf("      %-2s : %s\n", param_names[i], row)
    end

    # The Hessian is symmetric in exact arithmetic. This residual confirms the
    # 2nd-order perturbation formula is implemented correctly (both parameter
    # orderings of the cross term included); it should read ≈ 0 (machine precision).
    Hsub  = H_all[j1, active, active]
    asym  = maximum(abs.(Hsub .- Hsub')) / max(maximum(abs.(Hsub)), eps())
    @printf("\n   max relative asymmetry |H−Hᵀ|/|H| = %.2e   (→ 0 confirms symmetric Hessian)\n", asym)

    println("\n── Per-multiplier margins ─────────────────────────────────────────────────")
    println("    j      C₁           C₂          C₂/C₁")
    for j in eachindex(mu)
        r = isfinite(C1_vec[j]) && C1_vec[j] != 0 ? C2_vec[j]/C1_vec[j] : NaN
        @printf("   %2d   %9.4f    %9.4f    %7.3f\n", j, C1_vec[j], C2_vec[j], r)
    end

    println("\n── Worst-case robustness margin ───────────────────────────────────────────")
    @printf("   C₁ (first  order, linear)    = %.5f\n", C1)
    @printf("   C₂ (second order, quadratic) = %.5f\n", C2)
    @printf("   relative difference |C₂−C₁|/C₁ = %.2f %%\n", 100*abs(C2-C1)/abs(C1))

    ## ─── Monte-Carlo sanity check ──────────────────────────────────────────────
    # At margin C, scaling the uncertainty box by C should bring the worst-case
    # |μ| to ≈ 1. Sample the box corners + random points at scale C₁ and C₂ and
    # report the largest |μ| reached. A well-calibrated margin gives max|μ| ≈ 1.
    println("\n── Monte-Carlo verification (perturb by C·Δp, recompute |μ|) ───────────────")
    const N_random = 200
    active_v = (ζ_var, ω_var, kp_var, kd_var, τ_var)
    act_idx  = findall(!iszero, active_v)
    n_av     = length(act_idx)

    # History window must accommodate the largest perturbed τ.
    C_max = max(C1, C2)
    τ_max = τ_nom * (1 + C_max * τ_var)
    pm_mc = DDEProblem{true}(DelayPDOsc!, u0, h, (0.0, T_nom), p_point;
        constant_lags=[τ_max])
    KA_mc = (2, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=12, verbosity=0, maxiter=20))
    dp_mc = dynamic_problemSampled(pm_mc, Solver_args, τ_max;
        Historyresolution=25, zerofixpont=false, affineinteration=1, Krylov_arg=KA_mc)

    # Perturbed |μ| for a ratio vector r ∈ [-1,1]^5 scaled by C.
    function rho_at(r, C)
        rζ, rω, rkp, rkd, rτ = r
        p_pert = (ζ_nom + rζ*C*ζ_var*abs(ζ_nom),
                  ω_nom + rω*C*ω_var*abs(ω_nom),
                  kp_point + rkp*C*kp_var*abs(kp_point),
                  kd_point + rkd*C*kd_var*abs(kd_point),
                  τ_nom + rτ*C*τ_var*abs(τ_nom),
                  T_nom)
        return spectralradius(dp_mc; p=p_pert)
    end

    make_r(av) = (r = zeros(5); for (k,i) in enumerate(act_idx); r[i]=av[k]; end; ntuple(i->r[i],5))

    function scan_box(C)
        ρmax = 0.0
        # all corners
        for c in Iterators.product(ntuple(_->(-1.0,1.0), n_av)...)
            ρmax = max(ρmax, rho_at(make_r(c), C))
        end
        # random interior/edge samples
        for _ in 1:N_random
            ρmax = max(ρmax, rho_at(make_r(ntuple(_->2rand()-1, n_av)), C))
        end
        return ρmax
    end

    rho_at(make_r(ntuple(_->0.0,n_av)), 0.0)   # warmup
    ρmax_C1 = scan_box(C1)
    ρmax_C2 = scan_box(C2)
    @printf("   max|μ| over box scaled by C₁ = %.5f   (target ≈ 1.0)\n", ρmax_C1)
    @printf("   max|μ| over box scaled by C₂ = %.5f   (target ≈ 1.0)\n", ρmax_C2)
    println("   (closer to 1.0 = better-calibrated margin; <1 conservative, >1 optimistic)")
end

println("\n" * "="^78)
println(" Done.")
println("="^78)
