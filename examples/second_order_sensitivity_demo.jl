# Second-Order Sensitivity Analysis Demo
# Compares first-order (linear) vs second-order (quadratic) robustness margins
# for the delayed PD-controlled oscillator.
#
# Layout:
#   Panel 1: ρ = max|μ| spectral radius (stable region)
#   Panel 2: C₁ map (first-order linear)
#   Panel 3: C₂ map (second-order quadratic, vertex check)
#   Panel 4: |C₂ - C₁| / C₁  (relative difference)
#
# Validation:
#   - MDBM finds stability boundary ρ=1 (reference)
#   - MDBM finds C=C_validate boundary with C₁ (linear prediction)
#   - MDBM finds C=C_validate boundary with C₂ (quadratic prediction)
#   - MC MDBM: perturb params N times, overlay actual ρ=1 boundaries
#   - Count violations: how many MC boundaries cross each C=C_validate contour
#
# Theory: docs/second_order_sensitivity_theory.md

using Pkg; Pkg.activate(".")
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using GLMakie
using MDBM
using Statistics

# ## ─── Resolution parameters ────────────────────────────────────────────────────
# # HIGH-RESOLUTION (production)
# const RES_sweep_kp  = 160
# const RES_sweep_kd  = 161
# const RES_mdbm_seed = 12
# const RES_mdbm_iter = 5
# const RES_N_mc      = 25     # number of MC MDBM validations

 # Midlle-RESOLUTION (production)
 const RES_sweep_kp  = 40
 const RES_sweep_kd  = 41
 const RES_mdbm_seed = 7
 const RES_mdbm_iter = 4
 const RES_N_mc      = 10     # number of MC MDBM validations

# # LOW-RESOLUTION (fast preview — uncomment to switch)
#  const RES_sweep_kp  = 20
#  const RES_sweep_kd  = 21
#  const RES_mdbm_seed = 8
#  const RES_mdbm_iter = 3
#  const RES_N_mc      = 5

## ─── Model ────────────────────────────────────────────────────────────────────
function DelayPDOsc!(du, u, h, p, t)
    ζ, ω, kp, kd, τ, T = p
    hτ = h(p, t - τ)
    F  = 0.1 * cos(2π * t / T)
    du[1] = u[2]
    #Unstable with +ω^2, stable with -ω^2
    du[2] = +ω^2 * u[1] - 2ζ*ω * u[2] - kp * hτ[1] - kd * hτ[2] + F
end

const ζ_nom, ω_nom, kp_nom, kd_nom, τ_nom, T_nom = -0.02, 1.0, 2.0, 2.0, 0.5, 0.5
const p_init = (ζ_nom, ω_nom, kp_nom, kd_nom, τ_nom, T_nom)

u0 = @MArray [0.0, 0.0]
h(p, t) = @MArray [0.0, 0.0]
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-7)
probMapping = DDEProblem{true}(DelayPDOsc!, u0, h, (0.0, T_nom), p_init;
    constant_lags=[τ_nom])

include(joinpath(@__DIR__, "..", "src", "robust_control.jl"))

## ─── Uncertainty ───────────────────────────────────────────────────────────────
const ζ_var, ω_var, kp_var, kd_var, τ_var, T_var = 0.02, 0.00, 0.05, 0.05, 0.02, 0.0
const delta_p_nom = (ζ_var*abs(ζ_nom), ω_var*abs(ω_nom), kp_var*abs(kp_nom),
                     kd_var*abs(kd_nom), τ_var*abs(τ_nom), 0.0)
const C_validate = 1.0

## ─── Dynamic problems ─────────────────────────────────────────────────────────
Neig_spec = 2
KA_spec   = (Neig_spec, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig_spec+10, verbosity=0, maxiter=20))
dp_spec = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1, Krylov_arg=KA_spec)

# Shared Neig=12 dual dp for C₁ and C₂ sweeps
Neig_fast = 12
KA_fast   = (Neig_fast, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig_fast+20, verbosity=0, maxiter=20))
dp_dual = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1,
    Krylov_arg=KA_fast, perturbation_size=0.0)

## ─── Sweep grids ──────────────────────────────────────────────────────────────
const kpv = LinRange(0.5, 4.0, RES_sweep_kp)
const kdv = LinRange(-0.5, 4.0, RES_sweep_kd)

## ─── Sweep functions ──────────────────────────────────────────────────────────
function sweep_spec!(Smap)
    println("  Spectral radius sweep ($(length(kpv))×$(length(kdv)))...")
    Threads.@threads for j in 1:length(kpv)
        print(" ", j)
        for i in 1:length(kdv)
            Smap[i,j] = spectralradius(dp_spec;
                p=(ζ_nom, ω_nom, kpv[j], kdv[i], τ_nom, T_nom))
        end
    end
    println(" done.")
end

function sweep_C!(dp, Cmap, order)
    println("  C$order sweep ($(length(kpv))×$(length(kdv)))...")
    # NOTE: with_timeout is NOT used inside a Threads.@threads loop — async tasks
    # + Base.throwto across worker threads can deadlock. Direct call instead.
    Threads.@threads for j in 1:length(kpv)
        kp_v = kpv[j]; print(" ", j)
        for i in 1:length(kdv)
            kd_v = kdv[i]
            p_loc = (ζ_nom, ω_nom, kp_v, kd_v, τ_nom, T_nom)
            delta_loc = (ζ_var*abs(ζ_nom), ω_var*abs(ω_nom),
                         kp_var*abs(kp_v), kd_var*abs(kd_v),
                         τ_var*abs(τ_nom), 0.0)
            C_loc, _, mu_loc, _ = robustness_margin(dp, p_loc, delta_loc;
                verbose=false, sens_order=order)
            ρ = maximum(abs.(mu_loc))
            if isfinite(C_loc) && ρ < 1.0
                Cmap[i,j] = C_loc
            end
        end
    end
    println(" done.")
end

## ─── Compute maps ─────────────────────────────────────────────────────────────
Spec_map = zeros(length(kdv), length(kpv))
C1_map   = fill(NaN, length(kdv), length(kpv))
C2_map   = fill(NaN, length(kdv), length(kpv))

println("Warming up...")
spectralradius(dp_spec; p=p_init)
robustness_margin(dp_dual, p_init, delta_p_nom; verbose=false, sens_order=1)
robustness_margin(dp_dual, p_init, delta_p_nom; verbose=false, sens_order=2)

println("\n[1/4] Spectral radius map...")
t_spec = @elapsed sweep_spec!(Spec_map)

println("\n[2/4] First-order (C₁) map...")
t_C1 = @elapsed sweep_C!(dp_dual, C1_map, 1)

println("\n[3/4] Second-order (C₂) map...")
t_C2 = @elapsed sweep_C!(dp_dual, C2_map, 2)

# Apply authoritative stability mask
unstable = Spec_map .>= 1.0
C1_map[unstable] .= NaN
C2_map[unstable] .= NaN

println("\nTimings:")
println("  ρ sweep:  ", round(t_spec; digits=1), " s")
println("  C₁ sweep: ", round(t_C1;  digits=1), " s")
println("  C₂ sweep: ", round(t_C2;  digits=1), " s")
println("  Ratio C₂/C₁: ", round(t_C2/t_C1; digits=2), "×")

## ─── Plot ─────────────────────────────────────────────────────────────────────
fig = Figure(size=(2400, 750))
ax1 = GLMakie.Axis(fig[1,1], title="ρ = max|μ| ($(round(t_spec;digits=1)) s)", xlabel="kp", ylabel="kd")
ax2 = GLMakie.Axis(fig[1,2], title="C₁ — first order ($(round(t_C1;digits=1)) s)", xlabel="kp", ylabel="kd")
ax3 = GLMakie.Axis(fig[1,3], title="C₂ — second order ($(round(t_C2;digits=1)) s)", xlabel="kp", ylabel="kd")
ax4 = GLMakie.Axis(fig[1,4], title="|C₂−C₁|/C₁  (relative difference)", xlabel="kp", ylabel="kd")

# Spectral radius (stable region only)
Spec_s = copy(Spec_map); Spec_s[Spec_s .>= 1.0] .= NaN
hm1 = heatmap!(ax1, kpv, kdv, log.(Spec_s'), colormap=:inferno)
Colorbar(fig[2,1], hm1, vertical=false, label="log(ρ)")

# Shared color range from C₁ (the reference)
fv1 = filter(isfinite, C1_map)
cr  = isempty(fv1) ? (0.0,10.0) : (minimum(fv1), maximum(fv1))

hm2 = heatmap!(ax2, kpv, kdv, C1_map', colormap=:viridis, colorrange=cr)
Colorbar(fig[2,2], hm2, vertical=false, label="C₁")

hm3 = heatmap!(ax3, kpv, kdv, C2_map', colormap=:viridis, colorrange=cr)
Colorbar(fig[2,3], hm3, vertical=false, label="C₂")

# Relative difference
relDiff = abs.(C2_map .- C1_map) ./ max.(abs.(C1_map), eps.(C1_map))
relDiff[.!isfinite.(relDiff)] .= NaN
hm4 = heatmap!(ax4, kpv, kdv, relDiff', colormap=:plasma,
    colorrange=(0.0, min(maximum(filter(isfinite,relDiff)), 0.5)))
Colorbar(fig[2,4], hm4, vertical=false, label="|C₂−C₁|/C₁")

save("examples/second_order_sensitivity_demo_maps.png", fig)
display(fig)

## ─── MDBM: stability + C₁ + C₂ boundaries ───────────────────────────────────
kpv_mb = LinRange(0.5, 4.0, RES_mdbm_seed)
kdv_mb = LinRange(-.5, 4.0, RES_mdbm_seed)

# ── ρ=1 stability boundary ──────────────────────────────────────────────────
println("\nMDBM ρ=1 boundary...")
foo_rho(kp,kd) = log(spectralradius(dp_spec;
    p=(ζ_nom,ω_nom,kp,kd,τ_nom,T_nom)))
mdbm_rho = MDBM_Problem(foo_rho, [MDBM.Axis(kpv_mb,"kp"), MDBM.Axis(kdv_mb,"kd")])
t_rho = @elapsed begin
    MDBM.solve!(mdbm_rho, RES_mdbm_iter, interpolationorder=0, verbosity=0)
    MDBM.solve!(mdbm_rho, 0, interpolationorder=1, verbosity=0)
end
x_rho, y_rho = getinterpolatedsolution(mdbm_rho)
println("  ρ=1: $(length(x_rho)) pts in $(round(t_rho;digits=1)) s")

# ── C₁ = C_validate MDBM boundary ──────────────────────────────────────────
function C1_func(kp,kd)::Float64
    p_loc  = (ζ_nom, ω_nom, kp, kd, τ_nom, T_nom)
    dp_loc = (ζ_var*abs(ζ_nom), ω_var*abs(ω_nom),
              kp_var*abs(kp), kd_var*abs(kd), τ_var*abs(τ_nom), 0.0)
    mus, s_nom = affine(dp_dual; p=p_loc)
    all_mu = ComplexF64.(mus[1]); ne = min(Neig_fast, length(all_mu))
    mu_n = all_mu[1:ne]; ev = mus[2]
    s1 = compute_sensitivities_dual(dp_dual, p_loc, dp_loc, mu_n, s_nom, ev)
    dv = collect(Float64, dp_loc)
    C  = minimum(map(eachindex(mu_n)) do j
        d = sum(abs(s1[j,i])*dv[i] for i in eachindex(dv))
        d < eps() ? Inf : (1.0 - abs(mu_n[j])) / d
    end)
    return C - C_validate
end

# ── C₂ = C_validate MDBM boundary ──────────────────────────────────────────
function C2_func(kp,kd)::Float64
    p_loc  = (ζ_nom, ω_nom, kp, kd, τ_nom, T_nom)
    dp_loc = (ζ_var*abs(ζ_nom), ω_var*abs(ω_nom),
              kp_var*abs(kp), kd_var*abs(kd), τ_var*abs(τ_nom), 0.0)
    mus, s_nom = affine(dp_dual; p=p_loc)
    all_mu = ComplexF64.(mus[1]); ne = min(Neig_fast, length(all_mu))
    mu_n = all_mu[1:ne]; ev = mus[2]
    s1, h2, _ = compute_second_order_sens(dp_dual, p_loc, dp_loc, mu_n, s_nom, ev)
    dv = collect(Float64, dp_loc)
    C  = minimum(map(eachindex(mu_n)) do j
        _worst_case_C_2nd(abs(mu_n[j]), s1[j,:], h2[j,:,:], dv)
    end)
    return C - C_validate
end

println("MDBM C₁ = $C_validate boundary...")
C1_func(kp_nom, kd_nom)  # warmup
mdbm_C1 = MDBM_Problem(C1_func, [MDBM.Axis(kpv_mb,"kp"), MDBM.Axis(kdv_mb,"kd")])
t_C1_mdbm = @elapsed begin
    MDBM.solve!(mdbm_C1, RES_mdbm_iter, interpolationorder=0, verbosity=0)
    MDBM.solve!(mdbm_C1, 0, interpolationorder=1, verbosity=0)
end
x_C1m, y_C1m = getinterpolatedsolution(mdbm_C1)
println("  C₁=$C_validate: $(length(x_C1m)) pts in $(round(t_C1_mdbm;digits=1)) s")

println("MDBM C₂ = $C_validate boundary...")
C2_func(kp_nom, kd_nom)  # warmup
mdbm_C2 = MDBM_Problem(C2_func, [MDBM.Axis(kpv_mb,"kp"), MDBM.Axis(kdv_mb,"kd")])
t_C2_mdbm = @elapsed begin
    MDBM.solve!(mdbm_C2, RES_mdbm_iter, interpolationorder=0, verbosity=0)
    MDBM.solve!(mdbm_C2, 0, interpolationorder=1, verbosity=0)
end
x_C2m, y_C2m = getinterpolatedsolution(mdbm_C2)
println("  C₂=$C_validate: $(length(x_C2m)) pts in $(round(t_C2_mdbm;digits=1)) s")

println("\nMDBM timing:  ρ=1 $(round(t_rho;digits=1))s  |  C₁ $(round(t_C1_mdbm;digits=1))s  |  C₂ $(round(t_C2_mdbm;digits=1))s  |  ratio C₂/C₁ $(round(t_C2_mdbm/t_C1_mdbm;digits=2))×")

## ─── MC + corner validation: collect first, plot once ────────────────────────
# We collect ALL random and corner boundary points into two arrays, then
# plot with a single legend entry each. No intermediate save/display inside the
# loop — a single final render shows everything cleanly.
val_vars   = (ζ_var, ω_var, kp_var, kd_var, τ_var)
active_idx = findall(!iszero, val_vars)
n_av       = length(active_idx)

println("\n[4/4] Validation: $RES_N_mc random + 2^$n_av=$(2^n_av) corner boundaries...")

const τ_max_mc = τ_nom * (1 + C_validate * τ_var)
pm_mc = DDEProblem{true}(DelayPDOsc!, u0, h, (0.0, T_nom), p_init;
    constant_lags=[τ_max_mc])
dp_mc = dynamic_problemSampled(pm_mc, Solver_args, τ_max_mc;
    Historyresolution=20, zerofixpont=false, affineinteration=1, Krylov_arg=KA_spec)

function make_r_mc(av)
    r = zeros(5); for (k,i) in enumerate(active_idx); r[i]=av[k]; end
    return ntuple(i->r[i], 5)
end

function perturbed_boundary_mc(r)
    rζ,rω,rkp,rkd,rτ = r; sc = C_validate
    foo(kp,kd) = log(spectralradius(dp_mc; p=(
        ζ_nom+rζ*sc*ζ_var*abs(ζ_nom), ω_nom+rω*sc*ω_var*abs(ω_nom),
        kp+rkp*sc*kp_var*abs(kp), kd+rkd*sc*kd_var*abs(kd),
        τ_nom+rτ*sc*τ_var*abs(τ_nom), T_nom)))
    prob = MDBM_Problem(foo, [MDBM.Axis(kpv_mb,"kp"), MDBM.Axis(kdv_mb,"kd")])
    try
        MDBM.solve!(prob, RES_mdbm_iter, interpolationorder=0, verbosity=0)
        MDBM.solve!(prob, 0, interpolationorder=1, verbosity=0)
        return getinterpolatedsolution(prob)
    catch; return Float64[], Float64[]; end
end

function pt_inside_contour(kp_pts, kd_pts, C_map_2d, C_lev)
    count_in = 0
    for (kp,kd) in zip(kp_pts, kd_pts)
        j = clamp(round(Int, (kp - kpv[1]) / step(kpv)) + 1, 1, length(kpv))
        i = clamp(round(Int, (kd - kdv[1]) / step(kdv)) + 1, 1, length(kdv))
        c = C_map_2d[i,j]
        !isnan(c) && c > C_lev && (count_in += 1)
    end
    return count_in
end

# Wrapped in a function to avoid Julia's top-level soft-scope warnings.
function collect_random_boundaries()
    xs = Float64[]; ys = Float64[]; v1 = 0; v2 = 0; tot = 0
    for n in 1:RES_N_mc
        r = make_r_mc(ntuple(_ -> 2rand()-1, n_av))
        xb, yb = perturbed_boundary_mc(r)
        length(xb) == 0 && continue
        append!(xs, xb); append!(ys, yb)
        v1 += pt_inside_contour(xb, yb, C1_map, C_validate)
        v2 += pt_inside_contour(xb, yb, C2_map, C_validate)
        tot += length(xb)
        println("  random $n/$RES_N_mc: $(length(xb)) pts (total=$(length(xs)))")
    end
    return xs, ys, v1, v2, tot
end

function collect_corner_boundaries()
    xs = Float64[]; ys = Float64[]; v1 = 0; v2 = 0; tot = 0
    for (k, c) in enumerate(Iterators.product(ntuple(_ -> (-1.0,1.0), n_av)...))
        xb, yb = perturbed_boundary_mc(make_r_mc(c))
        length(xb) == 0 && continue
        append!(xs, xb); append!(ys, yb)
        v1 += pt_inside_contour(xb, yb, C1_map, C_validate)
        v2 += pt_inside_contour(xb, yb, C2_map, C_validate)
        tot += length(xb)
        println("  corner $k/$(2^n_av): $(length(xb)) pts")
    end
    return xs, ys, v1, v2, tot
end

mc_x, mc_y, violations_C1_mc, violations_C2_mc, total_mc = collect_random_boundaries()
corner_x, corner_y, violations_C1_cor, violations_C2_cor, total_cor = collect_corner_boundaries()

## ─── Final plot with all overlays and legends ────────────────────────────────
# MDBM contour lines
for ax in (ax1, ax2, ax3, ax4)
    scatter!(ax, x_rho, y_rho, color=:cyan,   markersize=3, label="ρ=1 (MDBM)")
end
scatter!(ax2, x_C1m, y_C1m, color=:red,  markersize=5, marker=:circle,  label="C₁=$(C_validate) (1st order, MDBM)")
scatter!(ax3, x_C2m, y_C2m, color=:lime, markersize=5, marker=:diamond, label="C₂=$(C_validate) (2nd order, MDBM)")
# Cross-overlay for comparison
scatter!(ax2, x_C2m, y_C2m, color=:lime, markersize=4, marker=:diamond, label="C₂=$(C_validate) (2nd order)")
scatter!(ax3, x_C1m, y_C1m, color=:red,  markersize=4, marker=:circle,  label="C₁=$(C_validate) (1st order)")

# C contour lines from heatmap (level = C_validate)
contour!(ax2, kpv, kdv, C1_map'; levels=[C_validate], color=:red,  linewidth=4, label="C₁=$(C_validate)")
contour!(ax3, kpv, kdv, C2_map'; levels=[C_validate], color=:blue, linewidth=3, label="C₂=$(C_validate)")

# MC random boundaries (all in one scatter call → one legend entry)
if length(mc_x) > 0
    for ax in (ax1, ax2, ax3)
        scatter!(ax, mc_x, mc_y, color=(:black, 0.65), markersize=4,
            label="MC random ($(RES_N_mc) runs)")
    end
end

# Corner boundaries (all in one scatter call → one legend entry)
if length(corner_x) > 0
    for ax in (ax1, ax2, ax3)
        scatter!(ax, corner_x, corner_y, color=(:orange, 0.65
        ), markersize=4,
            label="$(2^n_av) corners (worst-case)")
    end
end

# Legends
axislegend(ax1, position=:rt, framevisible=true)
axislegend(ax2, position=:rt, framevisible=true)
axislegend(ax3, position=:rt, framevisible=true)

save("examples/second_order_sensitivity_demo_maps.png", fig)
display(fig)

## ─── Summary ──────────────────────────────────────────────────────────────────
println("\n── Validation Summary ─────────────────────────────────────────────────────")
println("  Random MC ($RES_N_mc runs, $total_mc pts):")
println("    Violations C₁=$(C_validate): $violations_C1_mc  ($(round(100*violations_C1_mc/max(total_mc,1);digits=1))%)")
println("    Violations C₂=$(C_validate): $violations_C2_mc  ($(round(100*violations_C2_mc/max(total_mc,1);digits=1))%)")
println("  Corners ($(2^n_av) worst-case, $total_cor pts):")
println("    Violations C₁=$(C_validate): $violations_C1_cor  ($(round(100*violations_C1_cor/max(total_cor,1);digits=1))%)")
println("    Violations C₂=$(C_validate): $violations_C2_cor  ($(round(100*violations_C2_cor/max(total_cor,1);digits=1))%)")
println("\n  C₂ improvement: MC $(violations_C1_mc-violations_C2_mc)  |  corners $(violations_C1_cor-violations_C2_cor) fewer violations")
println("  If violations_C₂ < violations_C₁ → 2nd order gives a more accurate prediction.")
println("\nFigure saved: examples/second_order_sensitivity_demo_maps.png")


display(fig)