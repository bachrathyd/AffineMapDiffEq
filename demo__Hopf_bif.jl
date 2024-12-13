#Demo: Delayed Nonline Oscill with nonlinearity
5 + 5

using Revise
using DDE_mapping

using BenchmarkTools
using Plots
theme(:dark)#:vibrant:dracula:rose_pine
plotly()
#const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )
default(size=(500, 500), titlefont=(15, "times"), legendfontsize=13, guidefont=(12, :white), tickfont=(12, :orange), guide="x", framestyle=:zerolines, yminorgrid=true, fontfamily="Computer Modern", label=nothing, dpi=600)

#using Profile
using StaticArrays
using DifferentialEquations

using LinearAlgebra

# using MDBM

using FileIO
using Suppressor

# using MDBM

using FileIO
# Governing equation

function DelayedNonlineOscill(u, h, p, t)
    # Parameters
    ζ, δ, b, τ, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^3 + μ * u[2]^5
    #ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^2*sign(u[2])
    # Update the derivative vector
    SA[dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
ζ = 0.01         # damping coefficient
δ = 3.1#0.2          # nat. freq
b = -2.8#stable per.orbit
b = -3.15#unstable per.orbit
b = -1#Hopf starting...
b = 0.1#Hopf starting...
τ = 2pi#0.5#2pi          # Time delay
μ = 5.0
p = ζ, δ, b, τ, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation ---------------
#initial condition
u0 = SA[0.001, 0.0]
#history function
h(p, t) = SA[0.0; 0.0]

Tlongsim = 5000.2
Tend = 27.0
prob_long = DDEProblem{false}(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), p; constant_lags=[τ])
#SciMLBase.AutoSpecialize
#SciMLBase.NoSpecialize
#SciMLBase.FullSpecialize

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-7)#

@time sol = solve(remake(prob_long, p=(ζ, δ, b, τ, μ)); Solver_args...);
plot(sol)
##
#last period of the long simulation:
t_select_period = 0.0:0.01:Tend
t_select_delay = LinRange(0, τ, 200)
sol_period = sol(sol.t[end] .- t_select_period)
sol_delay = sol(sol.t[end] .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t, sol_period[1, :])
plot!(sol_period.t, sol_period[2, :])
plot!(sol_delay.t, sol_delay[1, :])
plot!(sol_delay.t, sol_delay[2, :])
#plot the phase space (u - du)
plot(sol_delay[1, :], sol_delay[2, :])
plot!(sol_period[1, :], sol_period[2, :])


## ---------------- simulation max amplitude ----------------------
## parameters

bv = [LinRange(-2.0, 1.95, 30)..., LinRange(1.99, 2.05, 30)...]
norm_solperiod = similar(bv)
@time Threads.@threads for ib in eachindex(bv)
    println(ib)
    bloc = bv[ib]
    #prob_long = DDEProblem(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), (ζ, δ, bloc, τ, μ); constant_lags=[τ])
    #sol = solve(prob_long; Solver_args...)#abstol,reltol

    sol = solve(remake(prob_long, p=(ζ, δ, bloc, τ, μ)); Solver_args...)#abstol,reltol

    #plot(sol)
    #last period of the long simulation:
    t_select_period = 0.0:0.01:Tend
    t_select_delay = eriod = 0.0:0.01:τ
    sol_period = sol(sol.t[end] .- t_select_period)
    sol_delay = sol(sol.t[end] .- t_select_delay)
    norm_solperiod[ib] = maximum(abs.(getindex.(sol_period.u, 1)))
end

plot(bv, norm_solperiod)
scatter!(bv, norm_solperiod)
##
#sol = solve(remake(prob_long, p=(ζ, δ, 2.010, τ, μ)); Solver_args...)#abstol,reltol
#plot(sol)





## ---------------- Affine mapping ---------------------------
#preparation for continuation:

using ForwardDiff
one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)

using KrylovKit
Neig = 12#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-25, krylovdim=3 + 20, verbosity=0, eager=true));

τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = Tend # timeperiod of the mapping

#Creating the problem
dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)

#bv_affine = LinRange(-2.0, 2.05, 201)#TODO: ez kell ha a balolbali szakaszt szeretném követeni
#bv_affine = LinRange(2.05, -2.0, 201)

bv_affine = [LinRange(-2.0, 1.95, 30)..., LinRange(1.99, 2.05, 30)...]
λ_μ₀ = Any[similar(bv_affine)...]

@time Threads.@threads for ib in eachindex(bv_affine)
    println(ib)
    bloc = bv_affine[ib]
    mu, saff, sol0 = affine(dp_0_Tfix; p=(ζ, δ, bloc, τ, μ))
    λ_μ₀[ib] = log.(abs.(mu[1])) / sol0.t[end]

end

plot(bv, norm_solperiod)
#scatter()
for k in 1:Neig
    plot!(bv_affine, getindex.(λ_μ₀, k))
end
plot!(legend=false)

## ---------------------------



##----------------------- DDE-Hopf amp stab -------------------


function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[2]
end
function affect_short!(integrator)
    if integrator.t > 0.1
        terminate!(integrator)
    end
end
cb_short = ContinuousCallback(condition, affect_short!, nothing)

Solver_args_T_short = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-5,
    :callback => cb_short)#

using ForwardDiff
#TODO: ez nem hiszem, hogy jó megoldás
Base.:convert(::Type{Float64}, x::ForwardDiff.Dual{Float64,Float64,1}) = x.value

using KrylovKit
Neig = 8#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-35, krylovdim=3 + 35, verbosity=0, eager=true));
#Creating the problem
Timeperiod = 20.0;#Tend



τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = Tend # timeperiod of the mapping

#Creating the problem
dp_0_cb = dynamic_problemSampled(prob_long, Solver_args_T_short, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=3,
    Krylov_arg=Krylov_arg)





#------------------ initalization in a critical point --------------------------
bcrit = bv_affine[findfirst(x -> x < 0, getindex.(λ_μ₀, 1))]
b_H_start = 0.1;#-1.0#bcrit - 0.05;

mu_c, saff_c, sol0_c = affine(dp_0_Tfix; p=(ζ, δ, bcrit, τ, μ));
mu_c[1]
T = mu_c[3]
As = mu_c[2]
Acrit = mu_c[2][1]




plot(sol_delay[1, :], sol_delay[2, :])
plot!(sol_period[1, :], sol_period[2, :])
plot!(getindex.(Acrit, 1), getindex.(Acrit, 2), lw=1)
#plot!(dp_0_Tfix.StateSmaplingTime, getindex.(Acrit,1),lw=2)
#plot!(dp_0_Tfix.StateSmaplingTime, getindex.(Acrit,2),lw=2)

Acrtit_cb = LinMap(dp_0_cb, Acrit; p=(ζ, δ, b_H_start, τ, μ))[1]
plot!(getindex.(Acrtit_cb, 1), getindex.(Acrtit_cb, 2), lw=2)





mu_c, saff_c, sol0_c = affine(dp_0_cb, Acrtit_cb * 10.0; p=(ζ, δ, b_H_start, τ, μ));
scatter!(getindex.(saff_c, 1), getindex.(saff_c, 2), lw=2)
scatter!(sol0_c[1, :], sol0_c[2, :])


# ----------- brute force naive continuation ----------------

#ustart = rand(typeof(dp_0_cb.Problem.u0), Nstep)
#ustart = ustart .* 0.0 .+ 1.5
ustart = saff_c
#bv_affine_H = LinRange(-2.0, -0.951233, 27)
#bv_affine_H = LinRange(-1.0, -4.0, 37)
#bv_affine_H = LinRange(1.2, 0.08, 27)
bv_affine_H = LinRange(b_H_start, -2.0, 150)
bv_affine_H = [LinRange(b_H_start, -1.8, 10)..., LinRange(-1.81, -1.96, 40)..., LinRange(-1.96, -1.98, 100)...]

bv_affine_H = [LinRange(b_H_start, 2, 100)...]
bv_affine_H = [b_H_start:0.01:2.0...]

bv_affine_H = [LinRange(b_H_start, 1.95, 30)..., LinRange(1.96, 2.01, 30)...]

λ_μ₀_Hopf = Any[similar(bv_affine_H)...]
Amp_H = Any[similar(bv_affine_H)...]
T_period_H = Any[similar(bv_affine_H)...]
#
#plot()
#Threads.@threads

@suppress_err begin
    @time for ib in eachindex(bv_affine_H)
        println(ib)
        @show bloc = deepcopy(bv_affine_H[ib])
        mu, saff, sol0 = affine(dp_0_cb, ustart; p=(ζ, δ, bloc, τ, μ))
        #mu, saff, sol0 = affine(dp_Hopf_callback; p=(ζ, δ, bloc, τ, μ))
        λ_μ₀_Hopf[ib] = log.(abs.(mu[1])) / sol0.t[end]
        #Amp[ib] = saff
        T_period_H[ib] = sol0.t[end]
        Amp_H[ib] = maximum(abs.(getindex.(saff, 1)))

        ustart = saff #TODO: enélkül nem áll be szépen (de csak azért nem csinálta, mert affineinteration nem érvényesült)

        #   plot!(sol0[1, :], sol0[2, :])
        #   aaa=plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=1, lab="")#
        #   display(aaa)


    end
end
#plot!(legend=false)

##--------
plot(bv, norm_solperiod)
marcolor = sign.(getindex.(λ_μ₀_Hopf, 1))
#, color=:bamako
#scatter!(bv_affine_H, Amp_H, zcolor=marcolor, lw=0, markerstrokewidth=0)
scatter!(bv_affine_H, Amp_H, lw=0, markerstrokewidth=0)
plot!(legend=false)
plot!(ylim=(-0.3, 1000.0))
plot!(ylim=())



for k in 1:Neig
    plot!(bv_affine_H, getindex.(λ_μ₀_Hopf, k))
    #for ib in eachindex(bv_affine_H)
    # scatter!([bv_affine_H[ib]], [λ_μ₀_Hopf[ib][k]])
    #end
end
plot!()
plot!(ylim=(-0.5, 1.0))

fig_bif = plot!(bv_affine_H, T_period_H ./ 10.0, lw=3)


# ------------------------------------------------------------------------------
# TODO: csak itt tartok -  tesztelgetés
fig_normerror = scatter(yaxis=:log)
fig_abs_mus = scatter(yaxis=:log)
fig_path = plot()
# Pseudo-archlength method
5 + 5
b_H_start = 0.2
λ_start = b_H_start
p = (ζ, δ, b_H_start, τ, μ)
p_start = p
#initialization of a step:


mu_c, s0_start, sol0_c = affine(dp_0_cb, Acrtit_cb * 10.0; p=p_start);
mu_c, s0_start, sol0_c = affine(dp_0_cb, s0_start; p=p_start);


Δλ = 0.05#05
#Δλ = 0.005
#Δλ = 0.05
p = p .+ (0.0, 0.0, Δλ, 0.0, 0.0)
mu_c, s0, sol0_c = affine(dp_0_cb, s0_start; p=p);
mu_c, s0, sol0_c = affine(dp_0_cb, s0; p=p);



scatter!(fig_bif, [p_start[3]], [maximum(abs.(getindex.(s0_start, 1)))], markersize=4, m=:cross)
scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=4, m=:cross)


#using JLD2
#@save "bifurcation_till_Fold_point.jld2" 
#@load "bifurcation_till_Fold_point.jld2"
norm_limit = 1e-3
iteration_gola=5
Niteration = iteration_gola
maxΔλ=0.2


@suppress_err begin
    λscale = 10.0
    lamscale = 1.0
    Nconti = 1
    for Nconti in 1:40#(5*340)#34
        println("---------------------------")
        @show Nconti
        Δu = s0 - s0_start
        @show Δλ = p[3] - λ_start

        s0_start = deepcopy(s0)
        p_start = deepcopy(p)
        @show λ_start = p[3]
        
        @show Niteration
        @show lamscale = 1.3 .^ (iteration_gola .- Niteration)
        lamscale=minimum([lamscale,1.5])

        @show lamscale=(Δλ * lamscale)>maxΔλ ? 1.0 : lamscale

        s0 .+= Δu * lamscale
        p = p .+ (0.0, 0.0, Δλ * lamscale, 0.0, 0.0)
        dp = dp_0_cb
        Nstep = size(dp.StateSmaplingTime, 1)

        scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=3)
        scatter!(fig_bif, xlim=(0.09, 0.35), ylim=(0.05, 0.15))
       
        pDual_direction=(0.0, 0.0, 1.0 , 0.0, 0.0)
        
        mus, s0, sol, p,Niteration,normerror=affine(dp, s0; p, pDual_dir=pDual_direction, Δu=Δu,Δλ_scaled=λscale*Δλ,norm_limit = 1e-3) ;

        eigval = mus[1]
        eigvec = mus[2]
        
        scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], m=:cross, markersize=4)
        scatter!(fig_abs_mus, p[3] .* ones(size(eigval)), abs.(eigval), markersize=3, m=:circle, yaxis=:log,xlabel="λ",ylabel="asb_mu")
        
        
        v0, sol = LinMap(dp, s0; p=p)
        plot!(fig_path,getindex.(v0, 1), getindex.(v0, 2))#,xlim=(),ylim=()
        plot!(fig_path,sol[1, :], sol[2, :],xlabel="u1",ylabel="u2",title="path", lw=2)

        scatter!(fig_normerror, [p[3]], [norm(s0 - v0)], markersize=5, m=:cross, yaxis=:log,xlabel="λ",ylabel="Norm error")
        scatter!(fig_bif, xlim=(0.00, 2.3), ylim=(0.00, 0.8),xlabel="λ",ylabel="norm(u)")
        out = plot(fig_bif, fig_normerror, fig_abs_mus,fig_path,layout=(2, 2),size=(900, 600))
        display(out)
    end
end

# using JLD2
# #@save "bifurcation_till_Fold_point.jld2" 
# @load "bifurcation_till_Fold_point.jld2"


scatter!(fig_bif, xlim=(0.00, 2.3), ylim=(0.00, 0.8))
scatter!(fig_normerror, xlim=(0.00, 2.3))
out = plot(fig_bif, fig_normerror, fig_abs_mus,fig_path,layout=(2, 2),size=(900, 600))
