function dynamic_problemSampled(prob, alg, maxdelay, Tperiod; Historyresolution=200, eigN=4, zerofixpont=true)
    StateSmaplingTime = LinRange(-maxdelay, Float64(0.0), Historyresolution)#TODO: Float64!!!
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
    #mus = eigsolve(s -> LinMap(dp, s; p=p), size(dp.StateSmaplingTime, 1), dp.eigN, :LM)
    Nstep = size(dp.StateSmaplingTime, 1)
    s_start = rand(typeof(dp.DDEdynProblem.u0), Nstep)

    #mus = eigsolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM,KrylovKit.Arnoldi(orth=KrylovKit.ClassicalGramSchmidt()))
   # mus = eigsolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM)
    # # vals, vecs, info = eigsolve(...) 

    mus =  getindex(schursolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM, KrylovKit.Arnoldi(krylovdim =dp.eigN*1+5,tol=1e-4,verbosity = 0)),[3,2,1])
#    mus =  getindex(schursolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM, orth=KrylovKit.ClassicalGramSchmidt()),[3,2,1])
 # mus =  getindex(schursolve(s -> LinMap(dp, s; p=p), s_start, dp.eigN, :LM,KrylovKit.Arnoldi()),[3,2,1])
    # T, vecs, vals, info = schursolve(...) with

    
    return mus[1], mus[2]#::Vector{ComplexF64}
end
function spectralradius(dp::dynamic_problemSampled; p=dp.DDEdynProblem.p)
    return Float64(abs(spectrum(dp; p=p)[1][1]))::Float64
end
function affine(dp::dynamic_problemSampled; p=dp.DDEdynProblem.p)
    #TODO: fixed dimension problem!!!!
    Nstep = size(dp.StateSmaplingTime, 1)
    s0 = zeros(typeof(dp.DDEdynProblem.u0), Nstep)
    #s0 = rand(typeof(dp.DDEdynProblem.u0), Nstep)#TODO - ez csak teszt
    v0 = LinMap(dp, s0; p=p)

    #println(norm(s0-v0))
    Nstep = size(dp.StateSmaplingTime, 1)
    s_start = rand(typeof(dp.DDEdynProblem.u0), Nstep)
    ####mus = eigsolve(s -> LinMap(dp, s + s0; p=p) - v0, s_start, dp.eigN, :LM)
    #####mus = eigsolve(s -> LinMap(dp, s + s0; p=p) - v0, size(s0, 1), dp.eigN, :LM)#, krylovdim=dp.eigN*2)
    
    #mus = getindex(schursolve(s -> LinMap(dp, s + s0; p=p) - v0, s_start, dp.eigN, :LM,orth::KrylovKit.ClassicalGramSchmidt()),[3,2,1])
    
    #mus = getindex(schursolve(s -> LinMap(dp, s + s0; p=p) - v0, s_start, dp.eigN, :LM, KrylovKit.Arnoldi()),[3,2,1])
    mus =  getindex(schursolve(s -> LinMap(dp, s+ s0; p=p) - v0, s_start, dp.eigN, :LM, KrylovKit.Arnoldi(krylovdim =dp.eigN*1+5,tol=1e-4,verbosity = 0)),[3,2,1])
    #TODO: schursolve

    s0 = real.(find_fix_pont(s0, v0, mus[1], mus[2]))


    ###println(norm(s0 - LinMap(dp, s0; p=p)))
    #TODO: it might be better to incluse the mus calcluations here too
    for k_fix_iteration in 1:10
        s0 = real.(find_fix_pont(s0, LinMap(dp, s0; p=p), mus[1], mus[2]))
        normerror=norm(s0 - LinMap(dp, s0; p=p))
             if (normerror)<1e-3 
                println("Norm of fixpont mapping: $normerror after : $k_fix_iteration itreation.")
             break
             end
    end

    #return mus[1]::Vector{ComplexF64}
    return mus, s0
end

function LinMap(dp::dynamic_problemSampled, s; p=dp.DDEdynProblem.p)# where T

    #s = [SA[sv[1+(k-1)*2], sv[2+(k-1)*2]] for k in 1:size(sv, 1)÷2]
    StateSmaplingTime = LinRange(-dp.maxdelay, Float64(0.0), size(s, 1))
    dt = StateSmaplingTime[2] - StateSmaplingTime[1]

    ###TODO: milyen interpoláció kell?
    ##itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
    ###itp = interpolate(s, BSpline(Linear()))
    ###itp = interpolate(s, BSpline(Linear()))
    ##Hist_interp_linear = scale(itp, StateSmaplingTime)
    ###    itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
    ###    Hist_inÖterp_linear = scale(itp, dp.StateSmaplingTime)
    ##hint(p, t) = Hist_interp_linear(t)

    hint(p, t) = interpolate_complex_on_grid(s, -dp.maxdelay, dt, t)
    #sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), h=hint,p=p), MethodOfSteps(BS3()))#, save_everystep=false)#abstol,reltol
    sol = solve(remake(dp.DDEdynProblem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg, adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol

#    sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(RK4()), adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol
   # sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(RK4()))#, save_everystep=false)#abstol,reltol
    #sol = solve(remake(dp.DDEdynProblem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(ABM43()), adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol
    #TODO: saveat=ts
    v = [getvalues(sol, ti) for ti in StateSmaplingTime .+ dp.Tperiod]
    #vv = reduce(vcat, v)
    #return vv::Vector{Float64}
    #return vv#::Vector{T}
    return v
end

#function LinMapPerturbed(dp::dynamic_problemSampled, sv::Vector{Float64})::Vector{Float64}
#    vv = LinMap(dp, s0 .+ sv) .- v0
#    return vv::Vector{Float64}
#end

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
function find_fix_pont(s0::AbstractVector, v0::AbstractVector, eigval, eigvec)
    x = (v0 - s0)

    ##AtA = conj( eigvec' .* eigvec)
    AtA = [eigvec[i]' * eigvec[j] for i in 1:size(eigvec, 1), j in 1:size(eigvec, 1)]
    Atx = [eigvec[i]' * x for i in 1:size(eigvec, 1)]
    vi = AtA \ Atx
    vi_mu = (vi .* ((eigval) ./ (eigval .- 1.0)))
    #A=transpose(mapreduce(permutedims, vcat, eigvec))
    #fix_v = v0 - A * (((A'A) \ (A' * x)) .* ((eigval) ./ (eigval .- 1.0)))

    fix_v = v0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, vi_mu))
    return fix_v
end

function interpolate_complex_on_grid(x0::AbstractVector, t0::Float64, dt::Float64, x::Float64)
    idx = clamp(floor(Int, (x - t0) / dt), 0, length(x0) - 2)
    t = (x - t0 - idx * dt) / dt
    return (1 - t) * x0[idx+1] + t * x0[idx+2]
end



