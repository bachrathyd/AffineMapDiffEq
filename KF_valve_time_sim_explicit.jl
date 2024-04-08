#Kádár Fanni - Szelep rezgés - sim_only
# Algeb. Delay Diff. Eq.
5 + 5


using Revise
#using DDE_mapping

using Interpolations

using LinearAlgebra
using BenchmarkTools
using Plots
plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LaTeXStrings
# using Memoization

using MDBM
#using MAT
using Peaks
using NaNStatistics

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

x = Float64[]
y = Float64[]
z = Float64[]


logdekmin = Float64[]
logdekminY1 = Float64[]
scatter(Float64[], Float64[], Float64[])
ζ = 0.39
δ = 3.0
τ = 0.25
#for β = LinRange(7.6,7.9,10)#7.82
for β = LinRange(7.5, 7.815, 30)#7.82
    #for β = LinRange(7.79, 7.815, 50)#7.82

    xx = Float64[]
    yy = Float64[]
    zz = Float64[]
    println("---------------")#7.82
    println(β)
    Φ = 48.2
    #q = 3.0;
    q = 10.0

    p_FIX = (ζ, δ, τ, β, Φ, q, "semmise")
    T = τ

    using NonlinearSolve

    u0_FIX = SA[2.0, 0.0, 2.5, q*Φ/2]
    valve_KF_DADE_exp_FIXPOINT(u0_FIX, p_FIX)
    probFIX = NonlinearSolve.NonlinearProblem(valve_KF_DADE_exp_FIXPOINT, u0_FIX, p_FIX)
    @show uFIX = NonlinearSolve.solve(probFIX, NewtonRaphson())



    for Nmutip in 1.55:-0.025:1.05#0.1:0.1:2
        println("Nmutip: $Nmutip")
        Base.:+(a::SVector, b::Bool) = a .+ b

        #h123(p, t) = SA[2.0, 0.0, 2.5]
        u0start = [4.42 * Nmutip, 0.0, 2.5 * Nmutip] .* 0.0 .+ uFIX[1:3] .* Nmutip#+rand(3).*Nmutip
        h123(p, t) = SA[u0start...]
        h4(p, t) = q * Φ / 2 * Nmutip * 0.0 + uFIX[4] * Nmutip
        #u0start = [0.0, Nmutip,  0.0] .+ uFIX[1:3]#+rand(3).*Nmutip
        #h123(p, t) = SA[u0start...]
        #h4(p, t) =uFIX[4]
        p = (ζ, δ, τ, β, Φ, q, h4)
        u0 = h123(p, 0.0)
        valve_KF_DADE_exp_fun = DDEFunction(valve_KF_DADE_exp)
        #kt=0
        #prob_valve_exp = DDEProblem(valve_KF_DADE_exp_fun, u0, h123, [0.0, T] .+ kt * τ, p, constant_lags=[τ])
        #integrator = init(prob_valve_exp, MethodOfSteps(Rodas5()),reltol=1e-8, abstol=1e-8)
        #
        #step!(integrator, τ) 
        #plot(integrator.sol)
        #
        #h44 = y4_fun(integrator.sol,p)
        #integrator.p = (ζ, δ, τ, β, Φ, q, h44)
        #check_error(integrator)

        t_final = [0.0]
        u_final = [u0]
        f_final = [h4(p, 0.0)]

        sizehint!(t_final, 1_000_000)
        sizehint!(u_final, 1_000_000)
        sizehint!(f_final, 1_000_000)
        kt = 0
        #kt=10000
        #for kt in 0:10000
        for _ in 0:1000
            kt += 1
            #if mod(kt, 50) == 0
            #    println(kt)
            #end
            prob_valve_exp = DDEProblem(valve_KF_DADE_exp_fun, u0, h123, [0.0, T] .+ kt * τ, p, constant_lags=[τ])

            #@time 
            sol_exp = solve(prob_valve_exp, MethodOfSteps(Rodas5()), reltol=1e-12, abstol=1e-12)
            #sol_exp = solve(prob_valve_exp, MethodOfSteps(BS3()),        reltol=1e-5, abstol=1e-5,dt=0.5)

            #plot!(sol_exp, linestyle=:dash)

            h44 = y4_fun(sol_exp, p)
            fser = [h44(p, t) for t in sol_exp.t]
            #PLOTFIG=plot!(sol_exp.t, fser, linestyle=:dash,xlim=(0,sol_exp.t[end]),linewidth = 2)

            p = (ζ, δ, τ, β, Φ, q, h44)
            u0 = sol_exp[end]

            append!(t_final, sol_exp.t)
            append!(u_final, sol_exp.u)
            append!(f_final, fser)

            if false#mod(kt, 100) == 0
                println(kt)
                y1 = getindex.(u_final, 1)
                #plot(t_final,y1,linewidth = 2)
                MaxInd = argmaxima(y1)
                y1locmax = y1[MaxInd]
                pp = scatter(t_final[MaxInd], y1[MaxInd], color=:red)
                display(pp)
            end
        end

        #file = matopen("matfile.mat", "w")
        #write(file, "u_final", u_final)
        #write(file, "f_final", f_final)
        #close(file)



        #plotly()
        # aa=       plot(t_final,getindex.(u_final,1),linewidth = 3)#, linestyle=:dash
        #        plot!(t_final,getindex.(u_final,2),linewidth = 3)#, linestyle=:dash
        #        plot!(t_final,getindex.(u_final,3),linewidth = 3)#, linestyle=:dash
        #        aa=plot!(t_final,f_final,linewidth = 3)#, linestyle=:dash
        #display(aa)
        Npoints = size(t_final, 1)
        t_final = t_final[floor(Int, Npoints / 5):end]
        u_final = u_final[floor(Int, Npoints / 5):end]
        f_final = f_final[floor(Int, Npoints / 5):end]

        y1 = getindex.(u_final, 1)
        y1_pert = y1 .- uFIX[1]

        #plot(t_final, y1_pert, linewidth=2)

        MaxInd = argmaxima(y1)
        y1locmax = y1_pert[MaxInd]


        logdek = log.(y1locmax[1:end-1] ./ y1locmax[2:end])
        #scatter!(t_final[MaxInd], y1[MaxInd])

        #scatter(t_final[MaxInd][2:end],log.(y1locmax[1:end-1]./y1locmax[2:end]),yticks=LinRange(0,1e-2,10))

        #2D
        #  PLOTaaa = scatter!(logdek,        y1locmax[1:end-1] .+ uFIX[1], ylabel="y1", xlabel="LogDek")
        #3D 
        #  PLOTaaa = scatter!(y1locmax[1:end-1]*0.0 .+ β, y1locmax[1:end-1] .+ uFIX[1],logdek,
        #      xlabel="β", ylabel="y1",zlabel="LogDek")
        append!(xx, y1locmax[1:end-1] * 0.0 .+ β)
        append!(yy, y1locmax[1:end-1] .+ uFIX[1])
        append!(zz, logdek)

        #xlim=(0, 1e-2))
        #ylim=(0, 7))

        #xticks=LinRange(0, 1e-2, 11),
        # display(PLOTaaa)

    end

    Imin = argmin(zz)
    Nsmmoth=5
    zzSmoth = movmean(zz, Nsmmoth)
    append!(logdekmin,zzSmoth[Imin])
    append!(logdekminY1,yy[Imin])

    append!(x, xx)
    append!(y, yy)
    append!(z, zz)
    #plot!([0,0.003],[1.0,1.0] .* uFIX[1], linestyle=:dash,linewidth = 3,title="β: $β",        xlim=(-0.01, 1e-2))

    theme(:dao)
    foo(x) = (x > 0) ? (x) + 0.02 : 0
    #foo(x)= x
    #  BBplot=scatter(x, y; zcolor=foo.(z))

    BBplot = scatter(y, z; zcolor=foo.(x))
    plot!= scatter(yy, zz)
    plot!(logdekminY1, logdekmin, linewidth=3)

    display(BBplot)
end


plotly()
BBplot = scatter(y, z; zcolor=foo.(x))
plot!(logdekminY1, logdekmin, linewidth=3)

plot()
