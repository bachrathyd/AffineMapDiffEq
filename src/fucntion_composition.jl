#(history)function compatibility for Krylovkit - incomplete

#julia> myfun(s1::String, s2::String) = invoke(myfun, tuple(Any, Any), s1, s2) # overwrite
#WARNING: Method definition myfun(String, String) [...] overwritten
#https://discourse.julialang.org/t/overwriting-functions/404



using LinearAlgebra
using BenchmarkTools
using StaticArrays


#h(p, t::Float64) = ...
#h(p, t::Float64, deriv::Type{Val{1}}) = ...
struct funcomp{T}#<: Function   F<:Function
    scaler::Vector{T}
    f::Vector{Any}#{Function}
    range::Vector{T}
end
function Base.show(io::IO, fc::funcomp) # I heve to remove <: Function to have pretty-print
    println(io, "scaler: $(fc.scaler)")
    println(io, "functinos: $(fc.f)")
    println(io, "range $(fc.range)")
end



# ------------- Initializations --------------------
function funcomp()
    funcomp{Float64}(Float64[], funcomp[], [0.0, 0.0])
end

function funcomp(foo::F) where {F<:Function}
    funcomp{Float64}([1.0], [foo], [0.0, 1.0])
end
function funcomp(foo::F, range::Vector{T}) where {F<:Function} where {T}
    funcomp{Float64}([1.0], [foo], range)
end


function funcomp(scaler::Float64, foo::F) where {F<:Function}
    funcomp{Float64}([scaler], [foo], [0.0, 1.0])
end
function funcomp(scaler::Float64, foo::F, range::Vector{T}) where {F<:Function} where {T}
    funcomp{Float64}([scaler], [foo], range)
end

function funcomp(scaler::Vector{Float64}, foos::Vector)# where {F<:Function}
    reduce!(funcomp{Float64}(scaler, Any[foos...], [0.0, 1.0]))
end

function funcomp(scaler::Vector{Float64}, foos::Vector{F}, range::Vector{T}) where {F<:Function} where {T}
    reduce!(funcomp{Float64}(scaler, Any[foos...], range))
end

function funcomp(foos::Vector)
    reduce!(funcomp{Float64}(ones(Float64, size(foos, 1)), Any[foos...], [0.0, 1.0]))
end
function funcomp(foos::Vector, range::Vector{T}) where {T}
    reduce!(funcomp{Float64}(ones(Float64, size(foos, 1)), Any[foos...], range))
end

function (fc::funcomp{A})(p::Any, t::Float64) where {A}# where {B}
    @inbounds mapreduce(x -> x[1] .* x[2](p, t), +, zip(fc.scaler, fc.f))
end


#--------------- Functions --------------------


#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::N, v::funcomp)::funcomp where {N<:Number}
    return funcomp(α .* v.scaler, v.f, v.range)
end
#foo = 10.0 * his
#foo([1,2], -1.3)

#function Base.:+(v::funcomp, w::funcomp)::funcomp
#    ????
#end


#Base.similar(v): a way to construct vectors which are exactly similar to v
function Base.similar(v::funcomp)::funcomp
    funcomp()
end
#TODO: this is  bug (or feature): foo([1,2], -1.3)

function empty!(w::funcomp)::funcomp
    deleteat!(w.f, eachindex(w.f))#TODO:it is might be not necessary, but it is safe
    deleteat!(w.scaler, eachindex(w.scaler))#TODO:it is might be not necessary, but it is safe
    #w.range .= 0.0
    return w
end



#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w
function LinearAlgebra.mul!(w::funcomp, v::funcomp, α::Number)::funcomp
    empty!(w)
    append!(w.f, v.f)
    append!(w.scaler, α .* v.scaler)
    w.range .= v.range
    return w
end


#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::funcomp, α::Number)::funcomp
    v.scaler .*= α
    return v::funcomp
end

#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::Number, v::funcomp, w::funcomp)::funcomp

    println("..............")
    println(α)
    
    println("-----------")
    println(size(w.f, 1))
    println(size(v.f, 1))
    println("-----------")
    println(size(w.scaler, 1))
    println(size(v.scaler, 1))
    println(size(w.range, 1))
    println(size(v.range, 1))


    append!(w.f, v.f)
    append!(w.scaler, α .* v.scaler)
    w.range[1] = maximum([v.range[1], w.range[1]])
    w.range[2] = minimum([v.range[2], w.range[2]])



    return reduce!(w)
end

#TODO: ezzel lenne jó convert(::Type{funcomp}, x::F) where F<:Function= funcomp(x)::funcomp
function LinearAlgebra.axpy!(α::Number, v::funcomp, fw::Function)::funcomp
    LinearAlgebra.axpy!(α, v, funcomp(fw, v.range))
end



#LinearAlgebra.axpby!(α, v, β, w): store in w the result of α*v + β*w
function LinearAlgebra.axpby!(α::Number, v::funcomp, β::Number, w::funcomp)::funcomp
    reduce!(axpy!(α, v, rmul!(w, β)))
end


function reduce!(fc::funcomp)::funcomp
    #println("-----------")
    #println(size(fc.f, 1))
    #println(size(fc.scaler, 1))
    fsunique = unique(fc.f)
    newscales = [sum(fc.scaler[fc.f .== fu]) for fu in fsunique]

    #println(size(fsunique, 1))
    #println(size(newscales, 1))
    empty!(fc)#range is not deleted!!!

    append!(fc.f, fsunique)
    append!(fc.scaler, newscales)

    #println(size(fc.f, 1))
    #println(size(fc.scaler, 1))
    return fc
end




#TODO: ehhez kell a struktúrába valami range (validity range-et definiálni)
using Integrals
#LinearAlgebra.dot(v,w): compute the inner product of two vectors
function LinearAlgebra.dot(v::funcomp, w::funcomp)::Float64
    #TODO: nem az egészet kell összeintegrálni, 
    #TODO: nem jók az idők, mert a függvény az csak a hsitory function. Csak a mappingben kell a sol és abból kivenni a cuccot...

    tstart = maximum([v.range[1], w.range[1]])
    tend = minimum([v.range[2], w.range[2]])

    prob = IntegralProblem((t, p) -> v(p, t)' * w(p, t), tstart, tend)
    VW = solve(prob, HCubatureJL()).u# reltol=1e-5, abstol=1e-5r
    #println(VW)
    return VW::Float64
end

#LinearAlgebra.norm(v): compute the 2-norm of a vector
function LinearAlgebra.norm(v::funcomp)::Float64
    return sqrt(dot(v, v))::Float64
end

#---------------- tests -------------------

foo1(p, t) = sin(t)
foo2(p, t) = 10.0
foo3(p, t) = t


his = funcomp(300.0, foo2, [-1.0, 10.0])
his = funcomp([300.0, 500], [foo1, foo2], [-1.0, 10.0])
foos = deepcopy(his)
empty!(foos)
foos
his


#TODO: this is  bug (or feature): foo([1,2], -1.3)
fooX = funcomp([5.0, 10.0], [foo3, sin], [-1.0, 0.5])
his = funcomp([300.0, 500], [foo1, foo2])
mul!(fooX, his, 10.0)


foo = similar(his)

his = funcomp([300.0, 500], [foo1, foo2], [1.0, 2.0])
rmul!(his, 0.01)
his
his([], -1.3)


his = funcomp([300.0, 500], [foo1, foo2], [-2.0, 0.5])
his2 = funcomp(foo3, [-1.0, 1.5])
axpy!(10.0, his, his2)
axpy!(10.0, his, foo2)



his = funcomp([300.0, 500], [foo1, foo2], [-2.0, 0.5])
his2 = funcomp(foo3, [-1.0, 1.5])
axpby!(10.0, his, 20.0, his2)
axpby!(10.0, his, 20.0, his2)


reduce!(his2)


his = funcomp([foo1, foo2], [-2.0, 0.5])
his1 = funcomp(foo3, [-0.0, 1.0])
his2 = funcomp(foo3, [-1.0, 1.5])

dot(his1, his2)

norm(his1)

#------- test DifEq. sol-al (is)-------------


his = funcomp([3.0, 0.1, 1.3], [foo1, foo2, foo3])
@benchmark his([1, 2], 1.5)
his = funcomp([3.0, 0.1, 1.0], [foo1, foo3, (p, t) -> sol(t)[2]])
@benchmark his([1, 2], 1.5)
scatter(1:0.1:125, [sol(t)[2] for t in 1:0.1:125])
plot!(t -> his([], t), range(0, 125, length=10000))


#------------------------ Mathieu test------------------


function f_now(p, t)
    ζ, δ, ϵ, b, τ, T = p
    # println("Computed $t")
    # println("Computed $p")
    SA[-(δ + ϵ * cos(t)), -2*ζ]
    #SA[-(δ+ϵ*sign(cos(t))), -2*ζ]
end
#@memoize 
function f_past(p, t)
    ζ, δ, ϵ, b, τ, T = p
    SA[b, 0]
end

function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ, ν, T = p
    dx = u[2]
    #ddx = [-(δ+ϵ*cos(t)), -2*ζ]' * u + b * h(p, t - τ)[1] # Traditional Delayed Mathieu equation
    #sign
    ddx = [-(δ + ϵ * (cos(t))), -2 * ζ]' * u + b * h(p, t - τ)[1]# + ν * h(p, t - τ, Val{1})[2]#The last pert is a test for neutral system
    SA[dx, ddx]
end


ζ = 0.03
δ = 1.1
ϵ = 0.2
τ = 2pi
b = 0.2
ν = 0.25 # Neutral coefficient
T = 10;
2pi;
p = ζ, δ, ϵ, b, τ, ν, T
#p = (ζ, ωn, k, τ,10.0)

u0 = SA[1.0, 1.0]

h(p, t::Float64) = SA[1.0; -0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T*1.0), p; constant_lags=[τ], neutral=true)

sol = solve(probMathieu, MethodOfSteps(BS3()))
plot(sol)

@benchmark solve(probMathieu, MethodOfSteps(BS3()))

Base.:+(a::SVector, b::Bool) = a .+ b
Nstep = 50
τmax = 2pi + 0.1
dpdp = dynamic_problemSampled(probMathieu, MethodOfSteps(BS3()), τmax, T; Historyresolution=Nstep, eigN=20, zerofixpont=true);

@time muaff, s0aff = affine(dpdp; p=p);
muaff[1]

@benchmark affine(dpdp; p=p)



#---------------------------Krylov-Functional-------------

his2 = funcomp(h, [-τmax, 0.0])
his2([], -1.0)

function LinMap(his::funcomp)::funcomp
    sol = solve(remake(probMathieu; u0=his(p, 0.0), tspan=(0.0, T), h=his, p=p), MethodOfSteps(BS3()); verbose=false)
    his = funcomp([1.0], [(p, t) -> sol(t + T)], his.range)
end

@benchmark hisout = LinMap(his2)

mus = getindex(schursolve(LinMap, his2, 5, :LM, KrylovKit.Arnoldi(krylovdim=18, tol=1e-90, verbosity=2)), [3, 2, 1])
mus = schursolve(LinMap, his2, 5, :LM)


## Different behavioure
T=0.3
T=6.0 # somehow the number of scalers and the number of function are not the same.
#TODO: question: what happens if the same function is used by differnet instance, and on is chaneg inplace. 
#especiall, it is racursive useage??!?! -> brute-forece solution: try deepcopy of all functions after it is provided.
T=10.0
