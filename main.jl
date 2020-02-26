# import Pkg; Pkg.add("JuMP")

include("jsr.jl")

A1 = [1 0; 1 0]
B1 = [0 1; 0 1]

A2 = [-1 0; 1 0]
B2 = [0 -1; 0 1]

v = Matrix{Float64}[]
push!(v, A2)
push!(v, B2)

lb, ub = jsr(v)

println("[", lb, ", ", ub, "]")


