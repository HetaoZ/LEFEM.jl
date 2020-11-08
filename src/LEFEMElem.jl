
module LEFEMElem
include("MathKits.jl")
using .MathKits
const mk = MathKits
using LinearAlgebra

# ---------------------------
# Constants
# ---------------------------
const CONSTRAIN_ALPHA = 1.0e20

# ---------------------------
# Classes
# ---------------------------
mutable struct Node
    id::Int
    x0::Vector{Float64}
    d::Vector{Float64}
    u::Vector{Float64}
    a::Vector{Float64}
end
Node(dim::Int) = Node(0, zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim), zeros(Float64, dim))
export Node

abstract type AbstractElem end
export AbstractElem

function create_elem_type(name)
    code = quote
        mutable struct $name <: AbstractElem
            id::Int
            elemtype::String
            link::Vector{Int} # link vector of node IDs       
            ptype::String # problem type: bar, beam, pstrain, pstress, sym, 3d
            dir::Int # sign of element volume, for failure checking
        end
        $name() = $name(-1,"none",Int[],"none",0)
        $name(elemtype, ptype) = $name(-1,elemtype,Int[],ptype,0)
    end
    return code
end

create_elem_type_and_eval(name) = eval(create_elem_type(name))


# ---------------------------
# Material definitions
# ---------------------------
include("material.jl")

# ---------------------------
# Element definitions
# ---------------------------
include("elem_tri3.jl")
export Tri3, elem_stress, elem_strain, elem_jacobi,  check_elem, integ_elem_brick


# ---------------------------
# Common Functions
# ---------------------------

function elem_x(elem::T, nodes::Array{Node}) where T <: AbstractElem
    dim = length(nodes[1].x0)
    n = length(elem.link)
    x = Matrix{Float64}(undef,n,dim)
    for i = 1:n, j = 1:dim
        x[i, j] = nodes[elem.link[i]].x0[j] + nodes[elem.link[i]].d[j]
    end
    return x
end

"""
Area of 2D elements
"""
function elem_area(elem,nodes)
    x = elem_x(elem,nodes)
    A = mk.polygon_area(x[:,1], x[:,2])
    return A
end

"""
Check if an element fails.
"""
function check_elem(elem)
    if sign(elem_jacobi(elem)) != sign(elem.dir)
        error("Element ( "*string(elem.id)*" ) failed.")
    end
end

###
end


