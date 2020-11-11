mutable struct Node
    id::Int
    x0::Vector{Float64}
    d::Vector{Float64}
    u::Vector{Float64}
    a::Vector{Float64}
end
Node(dim::Int) = Node(0, zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim))

mutable struct Convex
    id::Int # global ID, not the boundary's local ID
    link::Vector{Int}
    normal::Vector{Float64}
end
Convex() = Convex(0,Int[], Float64[])

