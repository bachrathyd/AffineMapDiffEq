#Kádár Fanni - Szelep rezgés - sim_only
# Algeb. Delay Diff. Eq.
5 + 5

using Revise
#using DDE_mapping

using Interpolations

using LinearAlgebra
using BenchmarkTools
using Plots
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings
# using Memoization

using MDBM
#using MAT
using Peaks
using NaNStatistics

using ProgressLogging

using LaTeXStrings

function valve_KF_DADE_exp(y, h, p, t)
    ζ, δ, τ, β, Φ, q, h4 = p
    y4 = -h4(p, t - τ) - (Φ^2.0 * y[1]^2.0) / 2.0 + Φ * y[1] * abs(2 * h4(p, t - τ) + y[3] + (Φ^2.0 * y[1]^2.0) / 4.0)^0.5#Singular "f=y4"

    d_y1 = y[2]
    d_y2 = -2 * ζ * y[2] - (y[1] + δ) + y[3] - y4 + h4(p, t - τ)
    d_y3 = β * (q - 1 / Φ * (y4 + h4(p, t - τ)))
    SA[d_y1, d_y2, d_y3]
end
function valve_KF_DADE_exp_FIXPOINT(y, p)
    ζ, δ, τ, β, Φ, q, h4 = p
    y4 = y[4]
    fX = -2 * y4 - (Φ^2.0 * y[1]^2.0) / 2.0 + Φ * y[1] * abs(2 * y4 + y[3] + (Φ^2.0 * y[1]^2.0) / 4.0)^0.5#Singular "f=y4"

    d_y1 = y[2]
    d_y2 = -2 * ζ * y[2] - (y[1] + δ) + y[3]
    d_y3 = β * (q - 1 / Φ * (y4 + y4))
    SA[d_y1, d_y2, d_y3, fX]
end



function y4_fun(sol, p)
    ζ, δ, τ, β, Φ, q, _ = deepcopy(sol.prob.p)
    y4series = [
        -sol.prob.p[end](p, t - τ) - (Φ^2.0 * sol[1, it]^2.0) / 2.0 +
        Φ * sol[1, it] * abs(2 * sol.prob.p[end](p, t - τ) + sol[3, it] + (Φ^2.0 * sol[1, it]^2.0) / 4.0)^0.5
        for (it, t) in enumerate(sol.t)
    ]
    #  return y4series
    nodes = (sol.t,)
    itp = extrapolate(interpolate(nodes, y4series, Gridded(Linear())), Flat())

    h4_interp(p, t) = itp(t)

    return h4_interp #, y4series
end
## - futtatás 

ζ = 0.39
δ = 3.0
τ = 0.25
β = 6
Φ = 48.2
#q = 3.0;
q = 10

#u0 = SA[0.1, 0.1, 0.1, 0.1,0.0]
#h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2,0.0]
h(p, t) = SA[2.0, 0.0, 2.5, q*Φ/2]


using NonlinearSolve

u0_FIX = SA[2.0, 0.0, 2.5, q*Φ/2]

p_FIX = (ζ, δ, τ, β, Φ, q, "semmise")
valve_KF_DADE_exp_FIXPOINT(u0_FIX, p_FIX)
probFIX = NonlinearSolve.NonlinearProblem(valve_KF_DADE_exp_FIXPOINT, u0_FIX, p_FIX)
@show uFIX = NonlinearSolve.solve(probFIX, NewtonRaphson())

Base.:+(a::SVector, b::Bool) = a .+ b



function run_MoS_sim(NNN,Nmul)
    T = τ
    h123(p, t) = SA[4.0*Nmul, 0.0, 2.5*Nmul]
    h4(p, t) = q * Φ / 2*Nmul

    p = (ζ, δ, τ, β, Φ, q, h4)
    u0 = h123(p, 0.0)
    valve_KF_DADE_exp_fun = DDEFunction(valve_KF_DADE_exp)


    t_final = [0.0]
    u_final = [u0]
    f_final = [h4(p, 0.0)]

    sizehint!(t_final, 1_000_000)
    sizehint!(u_final, 1_000_000)
    sizehint!(f_final, 1_000_000)
    kt = 0
    #kt=10000
    #for kt in 0:10000
    @progress  for kt in 0:NNN#2000
        #kt += 1

        prob_valve_exp = DDEProblem(valve_KF_DADE_exp_fun, u0, h123, [0.0, T] .+ kt * τ, p, constant_lags=[τ])

        #@time 
        # sol_exp = solve(prob_valve_exp, MethodOfSteps(Rodas5()), reltol=1e-12, abstol=1e-12)

        resol = 1e-12
        sol_exp = solve(prob_valve_exp, MethodOfSteps(BS3()), reltol=resol, abstol=resol, dt=0.5)

        # plot!(sol_exp, linestyle=:dash)

        h44 = y4_fun(sol_exp, p)
        fser = [h44(p, t) for t in sol_exp.t]

        # PLOTFIG=plot!(sol_exp.t, fser, linestyle=:dash,xlim=(0,sol_exp.t[end]),linewidth = 2)
        # display(PLOTFIG)
        p = (ζ, δ, τ, β, Φ, q, h44)
        u0 = sol_exp[end]

        append!(t_final, sol_exp.t)
        append!(u_final, sol_exp.u)
        append!(f_final, fser)

        if false #mod(kt, 100) == 0
            println(kt)
            y1 = getindex.(u_final, 1)
            #plot(t_final,y1,linewidth = 2)
            MaxInd = argmaxima(y1)
            y1locmax = y1[MaxInd]
            pp = scatter(t_final[MaxInd], y1[MaxInd], color=:red)
            display(pp)
        end
    end

    return (t_final, u_final, f_final)

end


#
ζ = 0.39
δ = 3.0
τ = 0.25
β = 6
Φ = 48.2
#q = 3.0;
q = 10

#u0 = SA[0.1, 0.1, 0.1, 0.1,0.0]
#h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2,0.0]

h(p, t) = SA[4.0, 0.0, 2.5, q*Φ/2]
@time (t_final, u_final, f_final) = run_MoS_sim(2000,1.5)
#plotly()
gr()

ASPR = 1.5

aa = plot(t_final, getindex.(u_final, 1), linecolor=:green, linewidth=2, size=(560, 280), label=L"y_1")
display(aa)
println("MoS-Done")

y1 = getindex.(u_final, 1)
y1_pert = y1 .- uFIX[1]

#plot(t_final, y1_pert, linewidth=2)

MaxInd = argmaxima(y1)
y1locmax = y1_pert[MaxInd]

foo(x) = x > 0 ? x : NaN

logdek = log.(foo.(y1locmax[1:end-1] ./ y1locmax[2:end]))
scatter!(t_final[MaxInd], y1[MaxInd], ylims=(1.5, 6.25),xlims=(0, maximum(t_final)), label=L"peaks", legend=false)

#scatter(t_final[MaxInd][2:end],log.(y1locmax[1:end-1]./y1locmax[2:end]),yticks=LinRange(0,1e-2,10))
scatter!([], [], color=:blue, label=L"Λ")
aa = scatter!(twinx(), t_final[MaxInd][2:end], logdek, color=:blue, ylims=[-0.001, 0.1],xlims=(0, maximum(t_final)), dpi=1000, legend=false)
display(aa)
savefig("KF_valve_sim_Fig2_A_beta_6.png")




## ------------------------------------------

ζ = 0.39
δ = 3.0
τ = 0.25
β =7.791
Φ = 48.2
#q = 3.0;
q = 10


#u0 = SA[0.1, 0.1, 0.1, 0.1,0.0]
#h(p, t) = SA[2.0,  0.0, 2.5,q*Φ/2,0.0]
h(p, t) = SA[4.0, 0.0, 2.5, q*Φ/2]
(t_final, u_final, f_final) = run_MoS_sim(20_000,1.5)
#plotly()
gr()

ASPR = 1.5
NDROP=50
aa = plot(t_final[1:NDROP:end], getindex.(u_final, 1)[1:NDROP:end], linecolor=:green, linewidth=2, size=(560, 280), label=L"y_1")
display(aa)
println("MoS-Done")

y1 = getindex.(u_final, 1)
y1_pert = y1 .- uFIX[1]

#plot(t_final, y1_pert, linewidth=2)

MaxInd = argmaxima(y1)
y1locmax = y1_pert[MaxInd]

foo(x) = x > 0 ? x : NaN



logdek = log.(foo.(y1locmax[1:end-1] ./ y1locmax[2:end]))
scatter!(t_final[MaxInd][1:NDROP:end], y1[MaxInd][1:NDROP:end], ylims=(1.5, 6.25), label=L"peaks",markerstrokewidth=0)

#scatter(t_final[MaxInd][2:end],log.(y1locmax[1:end-1]./y1locmax[2:end]),yticks=LinRange(0,1e-2,10))
scatter!([], [], color=:blue, label=L"Λ",markerstrokewidth=0)
aa = scatter!(twinx(), t_final[MaxInd][2:end][1:NDROP:end], logdek[1:NDROP:end], color=:blue, dpi=1000, legend=false,markerstrokewidth=0, ylims=[-0.0000, 0.0003])#
display(aa)

savefig("KF_valve_sim_Fig2_A_beta_7_791.png")

