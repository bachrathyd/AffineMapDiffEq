b = 0.2
p = ζ, δ, ϵ, b, τ, T

Nstepv =floor.(Int,exp.( LinRange(log(20), log(1e5), 21))) # Ez jól néz ki
Tmean = zeros(Float64, size(Nstepv, 1))
Tstd = zeros(Float64, size(Nstepv, 1))
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 4.0
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.10
for kk in 1:size(Nstepv, 1)
    runpercent=kk / size(Nstepv, 1)*100
    #Nstep = 100
    Nstep = Nstepv[kk]
    τmax = 2pi + 0.1

    print("  - run [%]: $runpercent, resolution [1]: $Nstep ,  ")

    #dpdp_Time = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=4, zerofixpont=true)

    dpdp_Time = dynamic_problemSampled(probMathieu, MethodOfSteps(RK4()), τmax, T; Historyresolution=Nstep, eigN=6, zerofixpont=true)

    t = @benchmark spectralradius($dpdp_Time; p=$p)
    print(t)
    Tmean[kk] = BenchmarkTools.median(t).time / 1e9
    Tstd[kk] = BenchmarkTools.std(t).time / 1e9
    println(" <<<<<<< ")
end

#plot(Nstepv,Tmean,grid=true,yerror=Tstd)
plot(Nstepv, Tmean, grid=false, ribbon=Tstd, fillalpha=0.5, label=false, yaxis=:log10, xaxis=:log10, color=:magenta)

for kC in (10.0 .^ (-12:-2))
    plot!(Nstepv,kC .* Nstepv, label=false)
    #plot!(Nstepv,kC .* Nstepv .^2 , label=false,linestyle=:dash)
end
plot!(Nstepv, Tmean, color=:blue, marker=(:circle, 3, 1.0), label="T_{CPU}",
xlabel="resolution [1]",ylabel="time [s]",yaxis=:log10, xaxis=:log10,
xticks = 10.0 .^ (-1:8),yticks = 10.0 .^ (-10:8),
xlims = (-maximum(-Nstepv)/3,maximum(Nstepv)*3),
ylims = (-maximum(-Tmean)/3,maximum(Tmean)*3),
grid=true,
gridlinewidth = 2)

