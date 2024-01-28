b = 0.2
p = ζ, δ, ϵ, b, τ, T
Nstepv = floor.(Int, exp.(LinRange(log(3), log(1e4), 20))) #21 # Ez jól néz ki
Tmean = zeros(Float64, size(Nstepv, 1))
Tstd = zeros(Float64, size(Nstepv, 1))
mumaxVec_error = zeros(Float64, size(Nstepv, 1))
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10#1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1#5.0
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.50#0.1

#the_alg=MethodOfSteps(ABM43())
the_alg = MethodOfSteps(BS3())
#the_alg=MethodOfSteps(RK4())
#the_alg=MethodOfSteps(Tsit5())

#dpdp_Time = dynamic_problemSampled(probMathieu, the_alg, τmax, T; Historyresolution=maximum(Nstepv) * 100, eigN=1, zerofixpont=true)
dpdp_Time = dynamic_problemSampled(probMathieu, the_alg, τmax, T; Historyresolution=maximum(Nstepv) * 10, eigN=1, zerofixpont=true)
@time mumaxFINE = spectralradius(dpdp_Time; p=p)
for kk in 1:size(Nstepv, 1)
    runpercent = kk / size(Nstepv, 1) * 100
    #Nstep = 100
    Nstep = Nstepv[kk]
    τmax = 2pi + 0.1


    dpdp_Time = dynamic_problemSampled(probMathieu, the_alg, τmax, T; Historyresolution=Nstep, eigN=1, zerofixpont=true)

    t = @benchmark spectralradius($dpdp_Time; p=$p)
    muMaxLoc = spectralradius(dpdp_Time; p=p)
    locerror=abs(muMaxLoc - mumaxFINE)
    mumaxVec_error[kk] = locerror
    Tmean[kk] = BenchmarkTools.median(t).time / 1e9
    Tstd[kk] = BenchmarkTools.std(t).time / 1e9

    
    locerror=mumaxVec_error[kk] 
    print("  - run [%]: $runpercent, resolution [1]: $Nstep ,  muerror: $locerror  ")
    print(t)

    println(" <<<<<<< ")
end

Xax2plot = Nstepv,
XLABEL="resolution [1]"
Xax2plot=mumaxVec_error
XLABEL="Mu error"
#plot(Nstepv,mumaxVec_error, yaxis=:log10, xaxis=:log10)


#plot(Xax2plot,Tmean,grid=true,yerror=Tstd)
plot!(Xax2plot, Tmean, grid=false, ribbon=Tstd, fillalpha=0.5,
    label=false, yaxis=:log10, xaxis=:log10, color=:gray)

for kC in (10.0 .^ (-12:5))
    plot!(Xax2plot, kC .* Xax2plot .^(-0.3333), label=false, yaxis=:log10, xaxis=:log10)
    #plot!(Nstepv, kC .* Nstepv .^(-2.0), label=false, yaxis=:log10, xaxis=:log10)
    #plot!(Xax2plot,kC .* Nstepv .^2 , label=false,linestyle=:dash)
end
plot!(Xax2plot, Tmean, color=:red, marker=(:circle, 3, 1.0), label="T_{CPU}",
    xlabel=XLABEL, ylabel="time [s]", yaxis=:log10, xaxis=:log10,
    xticks=10.0 .^ (-10:8), yticks=10.0 .^ (-10:8),
    xlims=(1e-10,1e1),#xlims=(-maximum(-Xax2plot) / 3, maximum(Xax2plot) * 3),
    ylims=(-maximum(-Tmean) / 3, maximum(Tmean) * 3),
    grid=true,
    gridlinewidth=2)

