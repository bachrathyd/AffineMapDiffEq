

#(history)function compatibility for Krylovkit - incomplete

#julia> myfun(s1::String, s2::String) = invoke(myfun, tuple(Any, Any), s1, s2) # overwrite
#WARNING: Method definition myfun(String, String) [...] overwritten
#https://discourse.julialang.org/t/overwriting-functions/404



using LinearAlgebra
using BenchmarkTools
#@benchmark  foo1([1,2], -1.3)
#@benchmark  his([1,2], -1.3)


#h(p, t::Float64) = ...
#h(p, t::Float64, deriv::Type{Val{1}}) = ...
struct funcomp{T}#<: Function   F<:Function
    scaler::Vector{T}
    f::Vector{Any}#{Function}
end
function Base.show(io::IO, fc::funcomp) # I heve to remove <: Function to have pretty-print
    print(io, "scaler: $(fc.scaler) , functinos: $(fc.f)")
end


function funcomp() 
    funcomp{Float64}(Float64[], funcomp[])
end
function funcomp(foo::F) where {F<:Function}
    funcomp{Float64}([1.0], [foo])
end
function funcomp(scaler::Float64, foo::F) where {F<:Function}
    funcomp{Float64}([scaler], [foo])
end
function funcomp(scaler::Float64, foo::F) where {F<:Function}
    funcomp{Float64}([scaler], [foo])
end

function funcomp(scaler::Vector{Float64}, foos::Vector{F}) where {F<:Function}
    reduce!(funcomp{Float64}(scaler, Any[foos...]))
end
function funcomp(foos::Vector)
    reduce!(funcomp{Float64}(ones(Float64,size(foos,1)), foos))
end

function (fc::funcomp{A})(p::Any, t::Float64) where {A}# where {B}
    @inbounds mapreduce(x -> x[1] .* x[2](p, t), +, zip(fc.scaler, fc.f))
end


foo1(p, t) = 10.0
foo2(p, t) = 10.0
foo3(p, t) = t

his = funcomp([300.0, 500,2.0], [foo1, foo2,foo1])
his2 = funcomp([300.0, 500], [foo1, his])
his3 = funcomp([foo1, foo2,foo1])

his = funcomp(10.0, foo3)
his = funcomp(foo3)
his((1.1, "c"), 2.0)



#p = (1.0, 2.0)
#h(p, t::Float64) = SA[10.0 * t; -10.0 * t]
#h(p, t::Float64, deriv::Type{Val{1}}) = SA[10.0, 10.0]


#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::N, v::funcomp)::funcomp where {N<:Number}
    return funcomp(α .* v.scaler, v.f)
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
foo = similar(his)
#TODO: this is  bug (or feature): foo([1,2], -1.3)

function empty!(w::funcomp)::funcomp
    deleteat!(w.f, eachindex(w.f))#TODO:it is might be not necessary, but it is safe
    deleteat!(w.scaler, eachindex(w.scaler))#TODO:it is might be not necessary, but it is safe
    return w
end

foos=deepcopy(his)
empty!(foos)
foos
his
#TODO: this is  bug (or feature): foo([1,2], -1.3)


#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w
function LinearAlgebra.mul!(w::funcomp, v::funcomp, α::Number)::funcomp
    empty!(w)
    append!(w.f, v.f)
    append!(w.scaler, α .* v.scaler)
    return w
end
fooX = funcomp([5.0], [foo3,sin])
his = funcomp([300.0, 500], [foo1, foo2])
mul!(fooX, his, 10.0)

#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::funcomp, α::Number)::funcomp
    v.scaler .*= α
    return v::funcomp
end
his = funcomp([300.0, 500], [foo1, foo2])
foo = rmul!(his, 0.01)
his([], -1.3)

#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::Number, v::funcomp, w::funcomp)::funcomp
    append!(w.f, v.f)
    append!(w.scaler, α .* v.scaler)
    return reduce!(w)
end
#TODO: ezzel lenne jó convert(::Type{funcomp}, x::F) where F<:Function= funcomp(x)::funcomp
function LinearAlgebra.axpy!(α::Number, v::funcomp, w::Function)::funcomp
    LinearAlgebra.axpy!(α, v, funcomp(w))
end
his = funcomp([300.0, 500], [foo1, foo2])
his2 = funcomp(foo3)
axpy!(10.0,his,his2)
axpy!(10.0,his,foo2)


#LinearAlgebra.axpby!(α, v, β, w): store in w the result of α*v + β*w
function LinearAlgebra.axpby!(α::Number, v::funcomp, β::Number, w::funcomp)::funcomp
    reduce!(axpy!(α,v,rmul!(w,β)))
end
his = funcomp([300.0, 500], [foo1, foo2])
his2 = funcomp(foo3)
axpby!(10.0,his,20.0,his2)
axpby!(10.0,his,20.0,his2)

function reduce!(fc::funcomp)::funcomp
    fsunique=unique(fc.f)
    newcales=[sum(fc.scaler[fc.f .== fu]) for fu in unique(fc.f)]
    empty!(fc)
    append!(fc.f, fsunique)
    append!(fc.scaler, newcales)
    return fc
end
reduce!(his2)



#TODO: ehhez kell a struktúrába valami range (validity range-et definiálni)
#using Integrals
##LinearAlgebra.dot(v,w): compute the inner product of two vectors
#function LinearAlgebra.dot(v::funcomp, w::funcomp)::Float64
#    #TODO: nem az egészet kell összeintegrálni, 
#    #TODO: nem jók az idők, mert a függvény az csak a hsitory function. Csak a mappingben kell a sol és abból kivenni a cuccot...
#    prob = IntegralProblem((t, p) -> v(p,t)' * w(p,t), v.t[1], v.t[end])
#    VW = solve(prob, HCubatureJL()).u# reltol=1e-5, abstol=1e-5r
#    return VW::Float64
#end
#
#push!(sol.u, sol.u[1])
##LinearAlgebra.norm(v): compute the 2-norm of a vector
#function LinearAlgebra.norm(v::ODESolution)::Float64
#    return sqrt(dot(v, v))::Float64
#end



his = funcomp([3.0, 0.1,1.3], [foo1, foo2,foo3])
@benchmark his([1,2],1.5)
his = funcomp([3.0, 0.1,1.0], [foo1, foo3,(p,t)->sol(t)[2]])
@benchmark his([1,2],1.5)
scatter([sol(t)[2] for t in 1:0.1:125])
plot(t->his([],t),range(0, 125, length=10000))
