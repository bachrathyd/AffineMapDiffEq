#ODEfunction compatibility for Krylovkit - incomplete

#Base.:*(α, v): multiply v with a scalar α, which can be of a different scalar type; in particular this method is used to create vectors similar to v but with a different type of underlying scalars.
function Base.:*(α::Number, v::ODESolution)
    vout = deepcopy(sol)
    vout.u .*= α
    vout.k .*= α
    return vout::ODESolution
end
#Base.similar(v): a way to construct vectors which are exactly similar to v
function Base.similar(v::ODESolution)::ODESolution
    vout = deepcopy(v)
    vout.u .*= 0.0
    vout.k .*= 0.0
    return vout::ODESolution
end
#LinearAlgebra.mul!(w, v, α): out of place scalar multiplication; multiply vector v with scalar α and store the result in w

#LinearAlgebra.rmul!(v, α): in-place scalar multiplication of v with α; in particular with α = false, v is the corresponding zero vector
function LinearAlgebra.rmul!(v::ODESolution, α::Number)::ODESolution
    v.u .*= α
    v.k .*= α
    return v::ODESolution
end
function LinearAlgebra.rmul!(v::ODESolution, α::Bool)::ODESolution
    println(α)
    v.u .*= α#0.0
    v.k .*= α#0.0
    return v::ODESolution
end
#LinearAlgebra.axpy!(α, v, w): store in w the result of α*v + w
function LinearAlgebra.axpy!(α::Number, v::ODESolution, w::ODESolution)::ODESolution
    error("Not working, yet!")
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
    error("Not working, yet!")
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
