#(history)function compatibility for Krylovkit - incomplete
#h(p, t::Float64) = ...
#h(p, t::Float64, deriv::Type{Val{1}}) = ...
struct mutablefuncton <: Function
    h::Function
end
mutablefuncton(p,t)=mutablefuncton.h(p,t)

hfoo=mutablefuncton(h)
hfoo(1.1,0.9)

p=(1.0,2.0)
h(p, t::Float64) = SA[10.0*t; -10.0*t]
h(p, t::Float64, deriv::Type{Val{1}}) = SA[10.0, 10.0]


#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::Number, v::Function)
    vout(p,t)=α .* v(p,t)
    return vout::Function
end
foo=5.0*h
foo(p,-1.3)


#Base.similar(v): a way to construct vectors which are exactly similar to v
function Base.similar(v::Function)::Function
    ##vout(p,t)=zero(typeof(v(p,t)))
    vout(p,t)=0.0 .* v(p,t)
    return vout::Function
end
foo=similar(h)
foo(p,-1.3)



#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w
function LinearAlgebra.mul!(w::Function, v::Function, α::Number)::Function
    w(p,t)=0.0 .* v(p,t)
    return w::Function
end


#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::Function, α::Number)::Function
    vstor=deepcopy(v)
    v=α .* vstor
    return v::Function
end
h=[]
h(p, t::Float64) = SA[t; -t]
foo=rmul!(h,100.0)
h(p,-1.3)




#julia> myfun(s1::String, s2::String) = invoke(myfun, tuple(Any, Any), s1, s2) # overwrite
#WARNING: Method definition myfun(String, String) [...] overwritten
#https://discourse.julialang.org/t/overwriting-functions/404


function LinearAlgebra.rmul!(v::ODESolution, α::Bool)::ODESolution
    println(α)
    v.u .*= α#0.0
    v.k .*= α#0.0
    return v::ODESolution
end
#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::Number, v::ODESolution, w::ODESolution)::ODESolution
    
    u0_axpy= α .* v.prob.u0 .+ w.prob.u0
    h_axpy(p, t) = α * v.prob.h(p, t) .+ w.prob.h(p, t)
    if sol1.prob.neutral
    h_axpy(p, t, deriv::Type{Val{1}}) =  α * v.prob.h(p, t,Val{1}) + w.prob.h(p, t,Val{1}) 
    end
    sol_axpy = solve(remake(sol1.prob;u0=u0_axpy,h=h_axpy), MethodOfSteps(BS3()))#TODO use the original method
    
    for _ in 1:size(w.u,1)
        pop!(w.u)
        pop!(w.t)
        pop!(w.k)
    end
    for k in 1:size(sol_axpy.u,1)
        push!(w.u,sol1.u[k])
        push!(w.t,sol1.t[k])
        push!(w.k,sol1.k[k])
    end
    return w::ODESolution
end

#LinearAlgebra.axpby!(α, v, β, w): store in w the result of α*v + β*w
function LinearAlgebra.axpy!(α::Number, v::ODESolution, β::Number,w::ODESolution)::ODESolution   
    h_axpy(p, t) = α .* v.prob.h(p, t) + β .* w.prob.h(p, t)
    u0_axpy= α .* v.prob.u0 + β .* w.prob.u0
    if sol1.prob.neutral
    h_axpy(p, t, deriv::Type{Val{1}}) =  α * v.prob.h(p, t,Val{1}) + w.prob.h(p, t,Val{1}) 
    end
    sol_axpy = solve(remake(sol1.prob;u0=u0_axpy,h=h_axpy), MethodOfSteps(BS3()))#TODO use the original method
    w=sol_axpy
    return w::ODESolution
end
using Integrals
#LinearAlgebra.dot(v,w): compute the inner product of two vectors
function LinearAlgebra.dot(v::ODESolution, w::ODESolution)::Float64
    #TODO: nem az egészet kell összeintegrálni, 
    #TODO: nem jók az idők, mert a függvény az csak a hsitory function. Csak a mappingben kell a sol és abból kivenni a cuccot...
    prob = IntegralProblem((t, p) -> v(t)' * w(t), v.t[1], v.t[end])
    VW = solve(prob, HCubatureJL()).u# reltol=1e-5, abstol=1e-5r
    return VW::Float64
end

push!(sol.u, sol.u[1])
#LinearAlgebra.norm(v): compute the 2-norm of a vector
function LinearAlgebra.norm(v::ODESolution)::Float64
    return sqrt(dot(v, v))::Float64
end
