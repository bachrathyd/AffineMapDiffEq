b = 0.2
p = ζ, δ, ϵ, b, τ, T

Nstepv = floor.((2:5:40) .^ 4)
Nstepv = vcat(floor.((2:1:10) .^ 3),floor.((20:10:100) .^ 3))
#Nstepv=floor.((2:2:8).^2)
Tmean = zeros(Float64, size(Nstepv, 1))
Tstd = zeros(Float64, size(Nstepv, 1))
BenchmarkTools.DEFAULT_PARAMETERS.samples = 1
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.0050
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.20
for kk in 1:size(Nstepv, 1)
    print("  - run % ")
    print(kk / size(Nstepv, 1))
    print("  ------ t[ms]: ")
    #Nstep = 100
    Nstep = Nstepv[kk]
    τmax = 2pi + 0.1
    println(Nstep)
    #dpdp_Time = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=4, zerofixpont=true)

    dpdp_Time = dynamic_problemSampled(probMathieu, MethodOfSteps(RK4()), τmax, T; Historyresolution=Nstep, eigN=1, zerofixpont=true)

    t = @benchmark spectralradius($dpdp_Time; p=$p)
    print(t)
    Tmean[kk] = BenchmarkTools.median(t).time / 1e6
    Tstd[kk] = BenchmarkTools.std(t).time / 1e6
    println(" <<<<<<< ")
end

#plot(Nstepv,Tmean,grid=true,yerror=Tstd)
plot(Nstepv, Tmean, grid=false, ribbon=Tstd, fillalpha=0.5, label=false, yaxis=:log10, xaxis=:log10)
plot!(Nstepv, Tmean, color=:blue, marker=(:circle, 8, 1.0), label="T_{CPU}", yaxis=:log10, xaxis=:log10)
#plot(Nstepv,Tmean)

