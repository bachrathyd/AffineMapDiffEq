function dynamic_problemSampled(prob, alg, maxdelay, Tperiod; Historyresolution=200,
    zerofixpont=true, affineinteration=1, Krylov_arg=())
    #dt=maxdelay / Historyresolution

    StateSmaplingTime = LinRange(-maxdelay, 0.0, Historyresolution)#TODO: Float64!!!
    #eigs = zeros(ComplexF64, eigN)
    #eigsA = Vector{Vector{ComplexF64}}(undef,eigN)
    #eigsA = [zeros(ComplexF64, Historyresolution) for _ in 1:eigN]
    #fixpont = Vector{typeof(prob.u0)}
    #{ComplexF64,Int64,Float64}
    dynamic_problemSampled(prob, alg, maxdelay, Tperiod, zerofixpont, affineinteration,
        StateSmaplingTime, Krylov_arg)
end
#function remake(dp::dynamic_problemSampled, kwargs...)
#    DifferentialEquations.remake(dp.Problem, kwargs...)
#end



###u0=SA[1.0,2.0]
###u0=[SA[10,20],SA[1.0,2.0],4.0,2]
##x=[1.0,2.0]
##srand=randsimilar(x,4)
function randsimilar(x::AbstractArray, N::Int)::Vector{typeof(x)}
    xrand = [randsimilar(xi, N) for xi in x]
end
function randsimilar(x::SVector, N::Int)::Vector{typeof(x)}
    xrand = rand(typeof(x), N)
end


function spectrum(dp::dynamic_problemSampled; p=dp.Problem.p)
    #mus = eigsolve(s -> LinMap(dp, s; p=p)[1], size(dp.StateSmaplingTime, 1), dp.eigN, :LM)
    Nstep = size(dp.StateSmaplingTime, 1)
    #s_start=[dp.Problem.u0 for _ in 1:Nstep] #TODO:fill!!!
    s_start = rand(typeof(dp.Problem.u0), Nstep)


    #randsimilar!(s_start)
    #EIGEN BASED
    # mus = eigsolve(s -> LinMap(dp, s; p=p)[1], s_start, dp.eigN, :LM)
    # # vals, vecs, info = eigsolve(...)

    #ISSI BASED
    #mus = issi_eigen(dp::dynamic_problemSampled,p=p)

    #SCHUR BASED
    #mus = getindex(schursolve(s -> LinMap(dp, s; p=p)[1], s_start, dp.eigN, :LM, KrylovKit.Arnoldi(krylovdim=dp.eigN +dp.KrylovExtraDim, tol=dp.KrylovTol, verbosity=0)), [3, 2, 1])

    #    mus = getindex(schursolve(s -> LinMap(dp, s; p=p)[1], s_start, dp.eigN, :LM, KrylovKit.Arnoldi(verbosity=3)), [3, 2, 1])
    mus = getindex(schursolve(s -> LinMap(dp, s; p=p)[1], s_start, dp.Krylov_arg...), [3, 2, 1])
    # T, vecs, vals, info = schursolve(...) with


    return mus[1], mus[2]#::Vector{ComplexF64}
end
function spectralradius(dp::dynamic_problemSampled; p=dp.Problem.p)
    if dp.zerofixpont
        return Float64(maximum(abs.(spectrum(dp; p=p)[1])))::Float64
    else
        return Float64(maximum(abs.(affine(dp; p=p)[1][1])))::Float64
    end
end




function partialpart(xSA)#::SVector)
    bb = [x.partials[1] for x in xSA]
    return SA[bb...]
    #return MVector(bb...);
end


function valuepart(xSA)#::SVector)
    bb = [x.value[1] for x in xSA]
    return SA[bb...]
    #return MVector(bb...);
end

function affine(dp::dynamic_problemSampled, s0::T; p=dp.Problem.p, pDual_dir=p .* 0.0, Δu=s0 .* 0.0, Δλ_scaled=1.0, norm_limit=1e-10) where {T}
    # # #maximum(abs.(getindex.(saff_c, 1)))
    # # #maximum(abs.(getindex.(saff_c_m1, 1)))
    # # #Δu_Δλ = (saff_c .- saff_c_m1) ./ 0.01
    one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)
    Niteration = 0
    #TODO: fixed dimension problem!!!!
    v0_dual = LinMap(dp, s0; p=p .+ one_espilon_Dual .* pDual_dir)[1]
    v0 = valuepart.(v0_dual)
    dv0dλ = partialpart.(v0_dual)


    Finished_itertaion = 0
    do_more_iteration = true
    while do_more_iteration


        #println(norm(s0-v0))
        Nstep = size(dp.StateSmaplingTime, 1)
        s_start = rand(typeof(dp.Problem.u0), Nstep)


        #if true#~DODOAU
        # println("Float perturbation")
        # EPSI_TODO_REMOVE=1e-3
        # s_start .*=  EPSI_TODO_REMOVE
        # TheMapping(s::T) = (LinMap(dp, s + s0; p=p)[1] - v0)::T
        #else
        #println("Dual perturbation - it seems to be faster! ;-)")
        #one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)
        TheMapping(s::T) = partialpart.(LinMap(dp, s * one_espilon_Dual + s0; p=p)[1] - v0)::T
        #end



        # s_start = rand(typeof(dp.Problem.u0), Nstep) * ForwardDiff.Dual(0.0, 1.0)

        #TODO: in case of very high tolerance it will use all the krylovdim!!!
        #mus = getindex(schursolve(TheMapping, s_start, dp.eigN, :LM, KrylovKit.Arnoldi(krylovdim=dp.eigN +dp.KrylovExtraDim, tol=dp.KrylovTol, verbosity=0)), [3, 2, 1])
        #mus = getindex(schursolve(TheMapping, s_start, dp.eigN, :LM, KrylovKit.Arnoldi(verbosity=2)), [3, 2, 1])
        mus = getindex(schursolve(TheMapping, s_start, dp.Krylov_arg...), [3, 2, 1])

        eigval = mus[1]
        eigvec = mus[2]
        #  mus = issi_eigen(dp::dynamic_problemSampled,p=p)
        #TODO: schursolve
        #println(size(mus[1],1))


        #    a0 = real.(find_fix_pont(s0, v0, mus[1], mus[2]))::T
        a0, Δλ_loc = find_fix_pont(s0, v0, eigval, eigvec, dv0dλ, Δu, Δλ_scaled)#::T
        p = p .+ pDual_dir .* Δλ_loc

        #println("Fix point calculation ---------- Start")
        #TODO: it might be better to incluse the mus calcluations here too
        for k_fix_iteration in 1:30  #TODO:use input parameters for this with default value

            Niteration += 1
            # println("find_fix_pont_start")
            s0 = a0
            v0 = LinMap(dp, s0; p=p)[1]
            a0, Δλ_loc = find_fix_pont(s0, v0, eigval, eigvec, dv0dλ, Δu, Δλ_scaled)#::T#TODO: kell a real?
            p = p .+ pDual_dir .* Δλ_loc
            # println("find_fix_pont_end")
            normerror = norm(s0 - v0)
            if (normerror) < norm_limit #TODO:use input parameters for this with default value
                #   println("Norm of fixpont mapping: $normerror after : $k_fix_iteration itreation.")
                #  println("Fix point calculation ---------- End")
                break
            end
            # # #             scatter!(fig_normerror, [p[3]], [norm(s0 - v0)], markersize=1, yaxis=:log)
            # # #             scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=1)
        end
        v0 = LinMap(dp, s0; p=p)[1]
        a0, Δλ_loc = (find_fix_pont(s0, v0, mus[1], mus[2], dv0dλ, Δu, Δλ_scaled))#::T#TODO: kell a real?
        p = p .+ pDual_dir .* Δλ_loc
        # println("find_fix_pont_end")
        # s0 = find_fix_pont(s0, LinMap(dp, s0; p=p), mus[1], mus[2])

        #return mus[1]::Vector{ComplexF64}
        s0 = a0
        v0, sol = LinMap(dp, s0; p=p)
        normerror = norm(s0 - v0)
        # # #     scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=2)
        Finished_itertaion += 1
        do_more_iteration = Finished_itertaion < dp.affineinteration
        if !do_more_iteration
            return mus, s0::T, sol, p, Niteration, normerror
        end
    end
end

function affine(dp::dynamic_problemSampled; p=dp.Problem.p)
    #TODO: fixed dimension problem!!!!


    #global NNN=0
    Nstep = size(dp.StateSmaplingTime, 1)
    s0 = rand(typeof(dp.Problem.u0), Nstep)
    if dp.zerofixpont
        #s0 = zeros(typeof(dp.Problem.u0), Nstep)
        s0 .*= 0.0
    end

    ##affine(dp, s0; p=p)
    #println("----------------------")
    #@show s0
    #@show typeof(s0)
    #println("------<<<<<<<<<<<<<<<<<<<--------")

    muSFix = affine(dp, s0; p=p)#First iteration

    #println(NNN)
    return muSFix

end

function LinMap(dp::dynamic_problemSampled, s::T; p=dp.Problem.p) where {T}#::T # where T


    #global NNN +=1
    #s = [SA[sv[1+(k-1)*2], sv[2+(k-1)*2]] for k in 1:size(sv, 1)÷2]
    #StateSmaplingTime = LinRange(-dp.maxdelay, Float64(0.0), size(s, 1))
    StateSmaplingTime = dp.StateSmaplingTime
    #dt = StateSmaplingTime[2] - StateSmaplingTime[1]

    #TODO: milyen interpoláció kell? #"ez és a solver" minimuma dominálja a rendet
    itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
    #itp = interpolate(s, BSpline(Linear()))
    Hist_interp_linear = scale(itp, StateSmaplingTime)
    #    itp = interpolate(s, BSpline(Cubic(Line(OnGrid()))))
    #    Hist_inÖterp_linear = scale(itp, dp.StateSmaplingTime)
    hint(p, t) = Hist_interp_linear(t) #TODO: ha úgyis fix a lépls, akkor ez nem is kell!!!
    hint(p, t, deriv::Type{Val{1}}) = Interpolations.gradient(Hist_interp_linear, t)[1]
    #hint(p, t) = itp(t) #TODO: akkor ez is elég!!!


    #NewTimePoints = StateSmaplingTime .+ dp.Tperiod
    #####TODO: ez miért lassabb
    ##saveNewPostions=NewTimePoints[NewTimePoints .>= 0.0]
    ##savePastPostions=NewTimePoints[NewTimePoints .< 0.0]
    ##sol = solve(remake(dp.Problem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg, adaptive=false, dt=dt,saveat=saveNewPostions)#, save_everystep=false)#abstol,reltol
    ###sol = solve(remake(dp.Problem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg, saveat=saveNewPostions)#, save_everystep=false)#abstol,reltol
    ##History_foo(x)=h(p,x)
    ##Solpast=History_foo.(savePastPostions)
    ##v=vcat(Solpast,sol.u)
    ##return v


    #hint(p, t) = interpolate_complex_on_grid(s, -dp.maxdelay, dt, t)
    #sol = solve(remake(dp.Problem; u0=hint(p, 0.0), h=hint,p=p), MethodOfSteps(BS3()))#, save_everystep=false)#abstol,reltol

    #sol = solve(remake(dp.Problem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg; verbose=false)#, save_everystep=false)#abstol,reltol

    #println("zzzzzzzzzzzzzzzzzzzzz")
    #@show typeof(dp.Tperiod)
    #@show dp.Tperiod
    #@show hint(p, 0.0)
    #@show typeof(hint(p, 0.0))
    sol = solve(remake(dp.Problem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p); dp.alg...)#, adaptive=dp.adaptive, dt=dt; verbose=false,reltol=1e-7)#, save_everystep=false)#abstol,reltol
    #@show sol.t[end]

    ####TODO: az u0- az eleve jön a h ból mint default paramater, de ha a múltat máshogy táromom, akkor lehet, hogy meg kellene tartani.
    #### - NEM jó, mert ha definiálv van az u0 a felhasználó által, akkor azt nem módosítja és nem lesz jó!!!
    ####  sol = solve(remake(dp.Problem; tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg, adaptive=false, dt=dt; verbose=false)#, save_everystep=false)#abstol,reltol

    #    sol = solve(remake(dp.Problem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg;verbose=false, abstol=1e-10,reltol=1e-10)#, save_everystep=false)#abstol,reltol

    #sol = solve(remake(dp.Problem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg, abstol=1e-5,reltol=1e-5)#, save_everystep=false)#abstol,reltol
    #sol = solve(remake(dp.Problem; u0=hint(p, Float64(0.0)), tspan=(Float64(0.0), dp.Tperiod), h=hint, p=p), dp.alg, reltol=1e-15, save_everystep=false)#, save_everystep=false)#abstol,reltol

    #sol = solve(remake(dp.Problem; u0=hint(p, 0.0), h=hint,p=p), MethodOfSteps(BS3()))#, save_everystep=false)#abstol,reltol
    #sol = solve(remake(dp.Problem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(RK4()), adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol
    #sol = solve(remake(dp.Problem; u0=hint(p, 0.0), tspan=(0.0, dp.Tperiod), h=hint, p=p), MethodOfSteps(ABM43()), adaptive=false, dt=dt)#, save_everystep=false)#abstol,reltol
    #TODO: saveat=ts - ez lassabbnak tűnik!

    NewTimePoints = StateSmaplingTime .+ sol.t[end]
    v = [getvalues(sol, ti) for ti in NewTimePoints]
    #vv = reduce(vcat, v)
    #return vv::Vector{Float64}
    #return vv#::Vector{T}
    return v, sol
end

#function LinMapPerturbed(dp::dynamic_problemSampled, sv::Vector{Float64})::Vector{Float64}
#    vv = LinMap(dp, s0 .+ sv)[1] .- v0
#    return vv::Vector{Float64}
#end

function getvalues(sol::ODESolution, t::T) where {T<:Real}
    if t < 0.0
        sol.prob.h(sol.prob.p, t)::typeof(sol.prob.u0)
    elseif t == 0.0
        sol.prob.u0::typeof(sol.prob.u0)
    else
        sol(t)::typeof(sol.prob.u0)
    end
end

#function find_fix_pont(s0::AbstractVector, v0::AbstractVector, eigval::AbstractVector,eigvec::Vector{<:AbstractVector})
function find_fix_pont(s0::T, v0::T, eigval, eigvec) where {T}
    x = (v0 - s0)

    ##AtA = conj( eigvec' .* eigvec)
    #   AtA = [eigvec[i]' * eigvec[j] for i in 1:size(eigvec, 1), j in 1:size(eigvec, 1)]
    #   Atx = [eigvec[i]' * x for i in 1:size(eigvec, 1)]
    AtA = [dot(eigvec[i], eigvec[j]) for i in 1:size(eigvec, 1), j in 1:size(eigvec, 1)]

    # println("------------------------------------")
    # println(AtA)
    Atx = [dot(eigvec[i], x) for i in 1:size(eigvec, 1)]
    ci = AtA \ Atx

    #ci =  Atx #TODO: ez ugyan azt adja Schur esetén!!!
    ci_mu = (ci .* ((eigval) ./ (eigval .- 1.0)))#TODO: Szabad ezt csinálni, a Schur-nál, nem a sajátértékkel kellenen skálázni... (vagy az pont kiesik valós függvényeknél???)
    #A=transpose(mapreduce(permutedims, vcat, eigvec))
    #fix_v = v0 - A * (((A'A) \ (A' * x)) .* ((eigval) ./ (eigval .- 1.0)))

    fix_v = v0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
    return fix_v, 0
end

function find_fix_pont(s0::T, v0::T, eigval, eigvec, dv0dλ, Δu, Δλ_scaled::Tlam) where {T,Tlam}
    x = (v0 - s0)

    Atx = [dot(eigvec[i], x) for i in 1:size(eigvec, 1)]

    Atdvdlam = [dot(eigvec[i], dv0dλ) for i in 1:size(eigvec, 1)]

    ΔuAt = [dot(Δu, eigvec[i]) for i in 1:size(eigvec, 1)]

    T_jac = vcat(hcat(diagm(eigval .- 1.0), Atdvdlam),
        hcat(ΔuAt', Δλ_scaled))

    ci_arch = T_jac \ vcat(-Atx, 0.0)

    ci = ci_arch[1:end-1]
    Δλ_loc = real.(ci_arch[end])
    ci_mu = ci .* (eigval)
    #A=transpose(mapreduce(permutedims, vcat, eigvec))
    #fix_v = v0 - A * (((A'A) \ (A' * x)) .* ((eigval) ./ (eigval .- 1.0)))

    fix_v = real.(v0 + mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu)))


    return fix_v::T, Δλ_loc::Tlam
end



#function find_fix_pont(s0::AbstractVector, v0::AbstractVector, eigval::AbstractVector,eigvec::Vector{<:AbstractVector})
function Constratint_find_fix_pont(s0::AbstractVector, v0::AbstractVector, eigval, eigvec)
    x = (v0 - s0)

    #### ##AtA = conj( eigvec' .* eigvec)
    #### #   AtA = [eigvec[i]' * eigvec[j] for i in 1:size(eigvec, 1), j in 1:size(eigvec, 1)]
    #### #   Atx = [eigvec[i]' * x for i in 1:size(eigvec, 1)]
    #### AtA = [dot(eigvec[i], eigvec[j]) for i in 1:size(eigvec, 1), j in 1:size(eigvec, 1)]
    #### 
    #### # println("------------------------------------")
    #### # println(AtA)
    Atx = [dot(eigvec[i], x) for i in 1:size(eigvec, 1)]
    #### ci = AtA \ Atx

    ci = Atx #TODO: ez ugyan azt adja Schur esetén!!!
    ci_mu = (ci .* ((eigval) ./ (eigval .- 1.0)))#TODO: Szabad ezt csinálni, a Schur-nál, nem a sajátértékkel kellenen skálázni... (vagy az pont kiesik valós függvényeknél???)
    #A=transpose(mapreduce(permutedims, vcat, eigvec))
    #fix_v = v0 - A * (((A'A) \ (A' * x)) .* ((eigval) ./ (eigval .- 1.0)))

    fix_v = v0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
    return fix_v
end


function interpolate_complex_on_grid(x0::AbstractVector, t0::Float64, dt::Float64, x::Float64)
    idx = clamp(floor(Int, (x - t0) / dt), 0, length(x0) - 2)
    t = (x - t0 - idx * dt) / dt
    return (1 - t) * x0[idx+1] + t * x0[idx+2]
end

function issi_eigen(dp::dynamic_problemSampled; p=dp.Problem.p)
    Nstep = size(dp.StateSmaplingTime, 1)

    s0 = zeros(typeof(dp.Problem.u0), Nstep)
    ## if dp.zerofixpont
    ##     v0=s0
    ## else
    v0 = LinMap(dp, s0; p=p)[1]
    ## end
    S = [rand(typeof(dp.Problem.u0), Nstep) for _ in 1:dp.eigN]
    H = zeros(Float64, dp.eigN, dp.eigN)#For Schur based calculation onlyS
    for _ in 1:12
        #V = [LinMap(dp, Si; p=p)[1] for Si in S]
        V = [LinMap(dp, Si + s0; p=p)[1] - v0 for Si in S] #TODO: ez nem jó, mert
        StS = [S[i]' * S[j] for i in 1:size(S, 1), j in 1:size(S, 1)]
        StV = [S[i]' * V[j] for i in 1:size(S, 1), j in 1:size(S, 1)]
        H = StS \ StV

        FShurr = schur(H)
        #H==FShurr.vectors * FShurr.Schur * FShurr.vectors'

        #S=FE.vectors*V;
        S = [sum(FShurr.vectors[:, i] .* V) for i in 1:dp.eigN]
        S .= S ./ norm.(S)
        #norm.(S)


    end
    Eigvals = eigvals(H)
    p = sortperm(Eigvals, by=abs, rev=true)

    return Eigvals[p], S[p]
end
