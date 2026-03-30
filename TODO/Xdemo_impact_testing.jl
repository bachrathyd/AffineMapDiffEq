

m = 1 # mass
k = 1 # stiffness
c = 0.001
μ = 50
k2 = 2000#
g = 0.17#0.25#gap
Acc = -0.2#-0.77#0.2          # nat. freq
#Acc = 0.1#0.2          # nat. freq
τ = sqrt(2) / 2.1#2pi/10
Amp=0.15
#plot()
#for Acc in -0.2:0.01:0.2 #0.1:0.02:0.6
#begin 
    Acc =-0.8
    τ = sqrt(2) / 2.1#2pi/10
    τ = sqrt(2) / 3000.0#2pi/10
    T = 2pi * 0.95
    Ω = 1.1#1.30065#  0.96#2pi / T
    timepause = [[0.0, 0.0]]#Dirac-Delata Impulse
    p = (m, k, c, μ, k2, g, Acc, τ, Amp, Ω, timepause)
    # prob_long = solve(remake(prob_long, p=(m, k, c, μ, k2, g, Acc, τ, Amp, Ωloc, timepause)); Solver_args...)#abstol,reltol




    #Tlongsim = 2000
    Tlongsim = 200
    #Tlongsim = 40000#Ezzel már szép a nagyítási függvény
    #prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p, neutral=true, callback=cb_wall_event, constant_lags=[τ])#
    prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p, callback=cb_wall_event)#
    #prob_long = DDEProblem(Diff_oscill, u0, h, (0.0, Tlongsim), p)#; constant_lags=[τ]

    #Parameters for the solver as a Dict (it is necessary to collect it for later use)
    Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-15, :dtmax => 1e-1)#; constrained = true
    #Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-15)#
    #Solver_args = Dict(:alg => MethodOfSteps(Rosenbrock23()), :verbose => false, :reltol => 1e-15)#
    @time sol = solve(prob_long; Solver_args...)
    bbb=plot(sol,xlim=(Tlongsim-3*T,Tlongsim),dense=false)
    #bbb=plot(sol,xlim=(0,Tlongsim))
    # plot!(sol,dense=false)
    plot!(sol(sol.t, Val{1}, idxs=2))
    # scatter!(sol.t[1:end-1],sol[1,1:end-1])
    # scatter!(sol.t[1:end-1],sol[2,1:end-1])
    # scatter!((sol.t[1:end-1]+sol.t[2:end])/2,(sol[2,1:end-1]-sol[2,2:end]) ./ (sol.t[1:end-1]-sol.t[2:end]))
    # 
    # plot!(sol.t,[sol(tloc,continuity=:left,idxs=2)-sol(tloc,continuity=:right,idxs=2) for tloc in sol.t])
    # 
    # 
    # 
    # scatter(sol[1,1:end-1],sol[2,1:end-1])

    aaa = plot(sol[1, (3*end) ÷ 4:end-1], sol[2, (3*end) ÷ 4:end-1])
    ccc=plot(aaa,bbb)
    display(ccc)
#end
display(ccc)