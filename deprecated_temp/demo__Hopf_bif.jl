#Demo: Delayed Nonline Oscill with nonlinearity
5 + 5


using Revise
using DDE_mapping

#using BenchmarkTools
using Plots
theme(:dark)#:vibrant:dracula:rose_pine
plotly()
#const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )
default(size=(900, 900), titlefont=(15, "times"), legendfontsize=13, guidefont=(12, :white), tickfont=(12, :orange), guide="x", framestyle=:zerolines, yminorgrid=true, fontfamily="Computer Modern", label=nothing, dpi=600)

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

normpower = Inf
#normpower = 2

function DelayedNonlineOscill(u, h, p, t)
    # Parameters
    ζ, δ, b, τ, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^3 + μ * u[2]^5
    #ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^2*sign(u[2])
    # Update the derivative vector
    @MArray [dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
ζ = 0.01         # damping coefficient
δ = 3.1#0.2          # nat. freq
b = -2.8#stable per.orbit
b = -3.15#unstable per.orbit
b = -1#Hopf starting...
b = 0.1#Hopf starting...
#Hopf starting...
b = 0.0354#Hopf starting...0.035272963112167
τ = 2pi#0.5#2pi          # Time delay
μ = 5.0
p = [ζ, δ, b, τ, μ]
#p = (ζ, ωn, k, τ,10.0)

# test simulation ---------------
#initial condition
u0 = @MArray [0.001, 0.0]
#history function
h(p, t) = @MArray [0.0; 0.0]

Tlongsim = 25000.2
Tend = 27.0
prob_long = DDEProblem{false}(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), p; constant_lags=[τ])
#SciMLBase.AutoSpecialize
#SciMLBase.NoSpecialize
#SciMLBase.FullSpecialize

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-7)#

@time sol = solve(remake(prob_long, p=[ζ, δ, b, τ, μ]); Solver_args...);
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

bv = [LinRange(-1.95, 0.0, 100)..., LinRange(0.03, 0.05, 100)..., LinRange(0.06, 2.009, 100)...]
norm_solperiod = similar(bv)
@time Threads.@threads for ib in eachindex(bv)
    println(ib)
    bloc = bv[ib]
    #prob_long = DDEProblem(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), (ζ, δ, bloc, τ, μ); constant_lags=[τ])
    #sol = solve(prob_long; Solver_args...)#abstol,reltol

    sol = solve(remake(prob_long, p=[ζ, δ, bloc, τ, μ]); Solver_args...)#abstol,reltol

    #plot(sol)
    #last period of the long simulation:
    t_select_period = 0.0:0.01:Tend
    t_select_delay = eriod = 0.0:0.01:τ
    sol_period = sol(sol.t[end] .- t_select_period)
    sol_delay = sol(sol.t[end] .- t_select_delay)
    norm_solperiod[ib] = norm(sol_period.u, normpower)
end

plot!(bv, norm_solperiod)
##
#sol = solve(remake(prob_long, p=[ζ, δ, 2.010, τ, μ]); Solver_args...)#abstol,reltol
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


bv_affine = [LinRange(-1.95, 1.95, 100)..., LinRange(1.99, 2.009, 100)...]
bv_affine = bv
λ_μ₀ = Any[similar(bv_affine)...]

@time Threads.@threads for ib in eachindex(bv_affine)
    println(ib)
    bloc = bv_affine[ib]
    mu, saff, sol0 = affine(dp_0_Tfix; p=[ζ, δ, bloc, τ, μ])
    λ_μ₀[ib] = log.(abs.(mu[1])) / sol0.t[end]

end

plot(bv, norm_solperiod, ylim=(-0.8, 1.0))
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
    if integrator.t > 0.02
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
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-20, krylovdim=22, verbosity=0, eager=true));
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



# ----------- brute force naive continuation ----------------

bv_affine_H = [LinRange(0.1, 1.95, 50)..., LinRange(1.96, 2.009, 30)...]

λ_μ₀_Hopf = Any[similar(bv_affine_H)...]
Amp_H = Any[similar(bv_affine_H)...]
T_period_H = Any[similar(bv_affine_H)...]

mu, ustart, sol0 = affine(dp_0_cb; p=[ζ, δ, bv_affine_H[1], τ, μ])
@suppress_err begin
    @time for ib in eachindex(bv_affine_H)
        println(ib)
        bloc = deepcopy(bv_affine_H[ib])
        mu, saff, sol0 = affine(dp_0_cb, ustart; p=[ζ, δ, bloc, τ, μ])
        #mu, saff, sol0 = affine(dp_Hopf_callback; p=[ζ, δ, bloc, τ, μ])
        λ_μ₀_Hopf[ib] = log.(abs.(mu[1])) / sol0.t[end]
        #Amp[ib] = saff
        T_period_H[ib] = sol0.t[end]
        #Amp_H[ib] = maximum(abs.(getindex.(saff, 2)))
        Amp_H[ib] = norm(saff, normpower)

        ustart = saff #TODO: enélkül nem áll be szépen (de csak azért nem csinálta, mert affineinteration nem érvényesült)

        #   plot!(sol0[1, :], sol0[2, :])
        #   aaa=plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=1, lab="")#
        #   display(aaa)


    end
end
#plot!(legend=false)
##--------
plot(bv, norm_solperiod, ylim=(-0.8, 1.0))
marcolor = sign.(getindex.(λ_μ₀_Hopf, 1))
#, color=:bamako
#scatter!(bv_affine_H, Amp_H, zcolor=marcolor, lw=0, markerstrokewidth=0)
scatter!(bv_affine_H, Amp_H, lw=0, markerstrokewidth=0)
plot!(legend=false)
plot!(ylim=())



for k in 1:Neig
    plot!(bv_affine_H, getindex.(λ_μ₀_Hopf, k))
    #for ib in eachindex(bv_affine_H)
    # scatter!([bv_affine_H[ib]], [λ_μ₀_Hopf[ib][k]])
    #end
end
plot!()
plot!(ylim=(-0.8, 8.0))

fig_bif = plot!(bv_affine_H, T_period_H ./ 10.0, lw=3)


# ------------------------------------------------------------------------------
#TODO: verbosity keyword, for displaying data
# TODO: csak itt tartok -  tesztelgetés
# már egész jól működik, de sok a hard-coded rész
# TODO:az indítást meg kellene csinlni megfelelő bifurkácóis pontból indítással
# TODO: és kapcsolódó két sajátvektort úgy kombinálni, hogy kilégítse a a periodikus pályához tartozó Poin-carré metszet-nek megadott függvényt.
# TODO: istabil Hopf pálya indítás teszt
# TODO: tesztelni, különféle diffegyenelteke, esztergálás autós,...tesztelni, peridódikus pálya bifukációját (késleltetett rezonancia görge), Duffin cucc egyenletet ... duffing-chain?!?
# TODO: tesztelni a verziószámokat
##

@warn "Itt tartok! - valahogy el kellene tudni indítani a Hopf pontból merőlegesen, de valamiért nem megy..."
#------------------ initalization in a critical point --------------------------
using NonlinearSolve
#Finding the Hopf points:
foo_Hopf(b_NL_test, p_dumy) = abs(affine(dp_0_Tfix; p=[ζ, δ, b_NL_test, τ, μ])[1][1][1]) - 1
prob = IntervalNonlinearProblem(foo_Hopf, [0.0, 0.2])
prob = NonlinearProblem(foo_Hopf, 0.2)
bHopf = solve(prob, SimpleBroyden(), abstol=1e-4, maxiters=100).u


#bHopf=0.03527

mu, saff, sol0 = affine(dp_0_Tfix; p=[ζ, δ, bHopf, τ, μ])
@show lams = log.((mu[1])) / sol0.t[end] # predicting Hopf point based on the comlex conjugater pair of critical eigenvalue


#The corresponding eigen value are correct and the lengts is proper [-τ,0], but the phase condition is not fullfilled
#but with a single mapping with the solution with the pahse condition simulete it as long as to create a proper phase
s0_shur = mu[2][1] * 0.03# scaling is necessary, because the it hase unit norm
#s0_shur=mu[2][1]*0.08# scaling is necessary, because the it hase unit norm
norm(s0_shur)
norm(s0_shur, normpower)
v0_shur, sol = LinMap(dp_0_cb, s0_shur; p=[ζ, δ, bHopf, τ, μ])


plot(dp_0_cb.StateSmaplingTime, getindex.(s0_shur, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(s0_shur, 2), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(v0_shur, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(v0_shur, 2), xlim=(-7, 5))

#now it is mapped to istself (approximately), so it is a good initial guess
norm(v0_shur .- LinMap(dp_0_cb, v0_shur; p=[ζ, δ, bHopf, τ, μ])[1]) / norm(v0_shur)
s0_shur = v0_shur
v0_shur, sol = LinMap(dp_0_cb, s0_shur; p=[ζ, δ, bHopf, τ, μ])


plot(dp_0_cb.StateSmaplingTime, getindex.(s0_shur, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(s0_shur, 2), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(v0_shur, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(v0_shur, 2), xlim=(-7, 5))

plot()
pDual_direction = (0.0, 0.0, 1.0, 0.0, 0.0)
a0_fix = s0_shur

plot!(dp_0_cb.StateSmaplingTime, getindex.(a0_fix, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(a0_fix, 2), xlim=(-7, 5))
norm(a0_fix, normpower)
mus, a0_fix, sol_fix, p_fix, Niteration, normerror = affine(dp_0_cb, 0.5 .* a0_fix; p=[ζ, δ, bHopf, τ, μ], pDual_dir=pDual_direction, Δu=s0_shur, Δλ_scaled=0.01, norm_limit=1e-3)

#mu, a0_fix, sol0 = affine(dp_0_cb,s0_shur; p=[ζ, δ, bHopf, τ, μ]) #nem lehet direktben használni, mert épp a szinuláris pontot néznénk
norm(a0_fix, normpower)

abs.(mus[1])

plot!(sol_fix)
plot!(dp_0_cb.StateSmaplingTime, getindex.(a0_fix, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(a0_fix, 2), xlim=(-7, 5))

a0_fix = LinMap(dp_0_cb, a0_fix; p=[ζ, δ, bHopf, τ, μ])[1]
norm(a0_fix, normpower)

#plot(sol_fix)
plot!(dp_0_cb.StateSmaplingTime, getindex.(a0_fix, 1), xlim=(-7, 5))
plot!(dp_0_cb.StateSmaplingTime, getindex.(a0_fix, 2), xlim=(-7, 5))


@show p_fix

scatter!(fig_bif, [bHopf], [0.0], markersize=6)
scatter!(fig_bif, [p_fix[3]], [norm(a0_fix, normpower)], markersize=6)#, m=:cross

for mmpow in 0.1:0.1:5.0
    #mus, a0_fix, sol_fix, p_fix, Niteration, normerror = affine(dp_0_cb, mmpow .* s0_shur; p=[ζ, δ, bHopf, τ, μ], pDual_dir=pDual_direction, Δu=s0_shur, Δλ_scaled=0.01, norm_limit=1e-3)
    mus, a0_fix, sol_fix, p_fix, Niteration, normerror = affine(dp_0_cb, mmpow .* s0_shur; p=[ζ, δ, bHopf, τ, μ])
    mu, saff, sol0 = affine(dp_0_cb, ustart; p=[ζ, δ, bloc, τ, μ])
    scatter!(fig_bif, [p_fix[3]], [norm(a0_fix, normpower)], markersize=2, m=:cross)#
end
scatter!(fig_bif)

p_start = (ζ, δ, bHopf, τ, μ)
λ_start = p_start[3]
s0_start = a0_fix .* 0.0
p = p_fix
s0 = a0_fix
# ~~~~~~~~~~~~~~~~~ Pseudo-archlength method ~~~~~~~~~~~~~~~~~~~~~~~

function branch_plot(branch, dp, λdirection)

    fig_normerror = scatter(yaxis=:log)
    fig_abs_mus = scatter(yaxis=:log)
    fig_path = plot()
global fig_bif
    #fig_bif = scatter([b[4][3] for b in branch], [norm(b[2], normpower) for b in branch], markersize=5, m=:cross)
    fig_bif = scatter!(fig_bif,[b[4][3] for b in branch], [norm(b[2], normpower) for b in branch], markersize=5, m=:cross)
    plot!(ylim=(-1.5, 2), xlim=(-2, 2.5))
    # mus, s0, sol=affine(dp, s0; p=p);
    for b in branch
        λ_loc = sum(b[4] .* λdirection)
        eigval = b[1][1]
        s0 = b[2]
        scatter!(fig_abs_mus, λ_loc .* ones(size(eigval)), abs.(eigval), markersize=3, m=:circle, yaxis=:log, xlabel="λ", ylabel="asb_mu")


        v0, sol = LinMap(dp, s0; p=b[4])
        plot!(fig_path, getindex.(s0, 1), getindex.(s0, 2))#,xlim=(),ylim=()
        plot!(fig_path, sol[1, :], sol[2, :], xlabel="u1", ylabel="u2", title="path", lw=2)

        scatter!(fig_normerror, [λ_loc], [norm(s0 - v0)], markersize=5, m=:cross, yaxis=:log, xlabel="λ", ylabel="Norm error")
        scatter!(xlabel="λ", ylabel="norm(u)")
        #scatter!(fig_bif, xlim=(0.00, 2.3), ylim=(0.00, 10.8))

    end
    out = plot(fig_bif, fig_normerror, fig_abs_mus, fig_path, layout=(2, 2), size=(900, 600))

    display(out)

end


function one_step_conti!(branch, dp, λdirection)
    println("---------------------------")
    @show Nconti
    Δu = branch[end][2] .- branch[end-1][2]
    #Δλ = branch[end][4]' * λdirection - branch[end-1][4]' * λdirection
    @show Δλ = sum(branch[end][4] .* λdirection) - sum(branch[end-1][4] .* λdirection)#valid for tuples too

    @show Niteration =branch[end][5]
    @show normerror =branch[end][6]
    @show lamscale = 1.3 .^ (iteration_goal .- Niteration)
    lamscale = minimum([lamscale, 1.5])
    @show lamscale = (Δλ * lamscale) > maxΔλ ? 1.0 : lamscale

    s0_new = branch[end][2] .+ Δu * lamscale
    pnew = branch[end][4] .+ λdirection * Δλ


    output = affine(dp, s0_new; p=pnew, pDual_dir=λdirection, Δu=Δu, Δλ_scaled=λscale * Δλ, norm_limit=1e-3)
    issuccessful = true
    if issuccessful
        push!(branch, output)
    end

end


# # # # b_H_start = 0.2#bHopf#0.2
# # # # b_H_start = bHopf#0.2
# # # # 
# # # # λdirection = [0.0, 0.0, 1.0, 0.0, 0.0]
# # # # Δλ = 0.03#05
# # # # branch = Any[]
# # # # push!(branch, affine(dp_0_cb; p=[ζ, δ, b_H_start, τ, μ]));
# # # # mu_start = branch[end][1]
# # # # s0_shur = mu_start[2][1] * 0.03# scaling is necessary, because the it hase unit norm
# # # # norm(s0_shur)
# # # # push!(branch, affine(dp_0_cb, s0_shur; p=branch[end][4] .+ λdirection * Δλ));
# # # # 
# # # # branch_plot(branch, dp_0_cb, λdirection)
# # # # 



@warn "ami ki van kommentelve azzal fut, de az újjal nem"

branch = Any[]

b_H_start = 0.2#bHopf#0.2
λ_start = b_H_start
p = (ζ, δ, b_H_start, τ, μ)
p_start = p
#initialization of a step:

@warn "Itt tartok! - valahogy el kellene tudni indítani a Hopf pontból merőlegesen, de valamiért nem megy..."
#------------------ initalization in a critical point --------------------------
using NonlinearSolve
#Finding the Hopf points:
foo_Hopf(b_NL_test,p_dumy)=abs(affine(dp_0_Tfix; p=(ζ, δ, b_NL_test, τ, μ))[1][1][1])-1
prob = IntervalNonlinearProblem(foo_Hopf,[0.0,0.2])
#prob = NonlinearProblem(foo_Hopf, 0.2)
bHopf = solve(prob,abstol=1e-6,maxiters=100).u

push!(branch, affine(dp_0_cb; p=p_start));
Δλ = 0.03

p = p .+ λdirection .* Δλ
s0_start=branch[1][2]
 push!(branch, affine(dp_0_cb, s0_start;p=p));
 

branch_plot(branch, dp_0_cb, λdirection)

##
#using JLD2
#@save "bifurcation_till_Fold_point.jld2" 
#@load "bifurcation_till_Fold_point.jld2"
norm_limit = 1e-3
iteration_goal = 4
Niteration = iteration_goal
maxΔλ = 0.2
λscale = 10.0
lamscale = 1.0
Nconti = 1

@suppress_err for Nconti in 1:50#(5*340)#34
    one_step_conti!(branch, dp_0_cb, λdirection)
    branch_plot(branch, dp_0_cb, λdirection)
end
branch_plot(branch, dp_0_cb, λdirection)



one_step_conti!(branch, dp_0_cb, λdirection)
branch_plot(branch, dp_0_cb, λdirection)





















 b_H_start = 0.2#bHopf#0.2
 λ_start = b_H_start
 p = (ζ, δ, b_H_start, τ, μ)
 p_start = p
 #initialization of a step:
 
 
 mu_c, s0_start, sol0_c = affine(dp_0_cb; p=p_start);
 
Δλ = 0.03#05
#Δλ = 0.005
#Δλ = 0.05 


 p = p .+ (0.0, 0.0, Δλ, 0.0, 0.0)
 mu_c, s0, sol0_c = affine(dp_0_cb, s0_start;p=p);


plot(dp_0_cb.StateSmaplingTime,getindex.(s0_shur, 1),xlim=(-7,5))
plot!(dp_0_cb.StateSmaplingTime,getindex.(s0_shur, 2),xlim=(-7,5))
plot!(dp_0_cb.StateSmaplingTime,getindex.(v0_shur, 1),xlim=(-7,5))
plot!(dp_0_cb.StateSmaplingTime,getindex.(v0_shur, 2),xlim=(-7,5))

scatter!(fig_bif, [p_start[3]], [norm(s0_start, normpower)], markersize=3, m=:cross)
scatter!(fig_bif, [p[3]], [norm(s0, normpower)], markersize=3, m=:cross)
plot!(ylim=(-1.5,1),xlim=(-2,2.2))
#plot!(xlim=(-0.01,0.06),ylim=(-0.02,0.1))





@suppress_err begin
    λscale = 10.0
    lamscale = 1.0
    Nconti = 1
    for Nconti in 1:80#(5*340)#34
        println("---------------------------")
        @show Nconti
        Δu = s0 - s0_start
        @show Δλ = p[3] - λ_start

        s0_start = deepcopy(s0)
        p_start = deepcopy(p)
        @show λ_start = p[3]

        @show Niteration
        @show lamscale = 1.3 .^ (iteration_goal .- Niteration)
        lamscale = minimum([lamscale, 1.5])

        @show lamscale = (Δλ * lamscale) > maxΔλ ? 1.0 : lamscale

        s0 .+= Δu * lamscale
        p = p .+ (0.0, 0.0, Δλ * lamscale, 0.0, 0.0)
        dp = dp_0_cb
        Nstep = size(dp.StateSmaplingTime, 1)

        scatter!(fig_bif, [p[3]], [norm(s0, normpower)], markersize=3)
        scatter!(fig_bif, xlim=(0.09, 0.35), ylim=(0.05, 0.15))
        plot!(ylim=(-0.5, 2), xlim=(-2, 2))
        plot!(xlim=(-0.01, 0.06), ylim=(-0.02, 0.1))

        pDual_direction = (0.0, 0.0, 1.0, 0.0, 0.0)

        mus, s0, sol, p, Niteration, normerror = affine(dp, s0; p=p, pDual_dir=pDual_direction, Δu=Δu, Δλ_scaled=λscale * Δλ, norm_limit=1e-3)


        # mus, s0, sol=affine(dp, s0; p=p);

        eigval = mus[1]
        eigvec = mus[2]

        scatter!(fig_bif, [p[3]], [norm(s0, normpower)], m=:cross, markersize=4)
        scatter!(fig_abs_mus, p[3] .* ones(size(eigval)), abs.(eigval), markersize=3, m=:circle, yaxis=:log, xlabel="λ", ylabel="asb_mu")


        v0, sol = LinMap(dp, s0; p=p)
        plot!(fig_path, getindex.(v0, 1), getindex.(v0, 2))#,xlim=(),ylim=()
        plot!(fig_path, sol[1, :], sol[2, :], xlabel="u1", ylabel="u2", title="path", lw=2)

        scatter!(fig_normerror, [p[3]], [norm(s0 - v0)], markersize=5, m=:cross, yaxis=:log, xlabel="λ", ylabel="Norm error")
        scatter!(fig_bif, xlim=(0.00, 2.3), ylim=(0.00, 10.8), xlabel="λ", ylabel="norm(u)")
        out = plot(fig_bif, fig_normerror, fig_abs_mus, fig_path, layout=(2, 2), size=(900, 600))
        display(out)
    end
end

# using JLD2
# #@save "bifurcation_till_Fold_point.jld2" 
# @load "bifurcation_till_Fold_point.jld2"


scatter!(fig_bif, xlim=(-2.00, 2.3), ylim=(0.00, 10.0))
scatter!(fig_normerror, xlim=(0.00, 2.3))
out = plot(fig_bif, fig_normerror, fig_abs_mus, fig_path, layout=(2, 2), size=(900, 600))
