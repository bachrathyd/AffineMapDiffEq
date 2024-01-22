function dynamic_problemSampled(prob, alg, maxdelay, Tperiod; Historyresolution=200, eigN=4, zerofixpont=true)
    StateSmaplingTime = LinRange(-maxdelay, 0.0, Historyresolution)
    eigs = zeros(ComplexF64, eigN)
    #eigsA = Vector{Vector{ComplexF64}}(undef,eigN)
    #eigsA = [zeros(ComplexF64, Historyresolution) for _ in 1:eigN]
    #fixpont = Vector{typeof(prob.u0)}
    #{ComplexF64,Int64,Float64}
    dynamic_problemSampled(prob, alg, maxdelay, Tperiod, StateSmaplingTime, eigN, eigs, zerofixpont)
end
#function remake(dp::dynamic_problemSampled, kwargs...)
#    DifferentialEquations.remake(dp.DDEdynProblem, kwargs...)
#end

function spectrum(dp::dynamic_problemSampled; p=dp.DDEdynProblem.p)
    mus = eigsolve(s -> LinMap(dp, s; p=p), size(dp.StateSmaplingTime, 1), dp.eigN, :LM)
    return mus[1]::Vector{ComplexF64}
end
function affine(dp::dynamic_problemSampled; p=dp.DDEdynProblem.p)
    #TODO: fixed dimension problem!!!!
    Nstep=size(dp.StateSmaplingTime, 1)*2;
    s0 = zeros(Nstep,)
    v0 = LinMap(dpdp, s0; p=p)
    #println(norm(s0-v0))
    mus = eigsolve(s -> LinMap(dpdp, s + s0; p=p) - v0, size(s0, 1), dp.eigN, :LM)#, krylovdim=dp.eigN*2)
    
    s0 = real.(find_fix_pont(s0, v0, mus[1], mus[2]))
    ###println(norm(s0 - LinMap(dpdp, s0; p=p)))
    #TODO: it might be better to incluse the mus calcluations here too
    for _ in 1:5
        s0 = real.(find_fix_pont(s0, LinMap(dpdp, s0; p=p), mus[1], mus[2]))
       # println(norm(s0 - LinMap(dpdp, s0; p=p)))
    end

    #return mus[1]::Vector{ComplexF64}
    return mus,s0
end

function LinMap(dp::dynamic_problemSampled, sv; p=dp.DDEdynProblem.p)# where T

    s = [SA[sv[1+(k-1)*2], sv[2+(k-1)*2]] for k in 1:size(sv, 1)÷2]
    StateSmaplingTime = LinRange(-dp.maxdelay, 0.0, size(s, 1))
    dt = StateSmaplingTime[2] - StateSmaplingTime[1]

 #TODO: milyen interpoláció kell?
    itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
    #itp = interpolate(s, BSpline(Linear()))
    #itp = interpolate(s, BSpline(Linear()))
    Hist_interp_linear = scale(itp, StateSmaplingTime)
    #    itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
    #    Hist_interp_linear = scale(itp, dp.StateSmaplingTime)
    hint(p, t) = Hist_interp_linear(t)

    #hint(p, t) = interpolate_complex_on_grid(s, -dp.maxdelay, dt, t)

    #sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), h=hint,p=p), MethodOfSteps(BS3()))#, save_everystep=false)#abstol,reltol
    #sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(RK4()), adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol
    sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(ABM43()), adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol
    #TODO: saveat=ts
    v = [getvalues(sol, ti) for ti in StateSmaplingTime .+ dp.Tperiod]
    vv = reduce(vcat, v)
    #return vv::Vector{Float64}
    return vv#::Vector{T}
end

function LinMapPerturbed(dp::dynamic_problemSampled, sv::Vector{Float64})::Vector{Float64}
    vv = LinMap(dp, s0 .+ sv) .- v0
    return vv::Vector{Float64}
end

function getvalues(sol::ODESolution, t::Real)
    if t < 0.0
        sol.prob.h(sol.prob.p, t)
    elseif t == 0.0
        sol.prob.u0
    else
        sol(t)
    end
end



#function find_fix_pont(s0::AbstractVector, v0::AbstractVector, eigval::AbstractVector,eigvec::Vector{<:AbstractVector})
function find_fix_pont(s0::AbstractVector, v0::AbstractVector, eigval,eigvec)
    x=(v0-s0)
    A=transpose(mapreduce(permutedims, vcat, eigvec))
    fix_v=v0 - A *( ((A'A) \ ( A' * x)) .* ((eigval) ./ (eigval .- 1.0)))
end

function interpolate_complex_on_grid(x0::AbstractVector, t0::Float64, dt::Float64, x::Float64)
    idx = clamp(floor(Int, (x - t0) / dt), 0, length(x0) - 2)
    t = (x - t0 - idx * dt) / dt
    return (1 - t) * x0[idx+1] + t * x0[idx+2]
end



