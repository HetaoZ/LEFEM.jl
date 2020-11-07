"""
《有限元法》P282：H-W变分原理不知道如何实现。
暂时先用最小位能原理。
使用对角元素乘大数法实现位移约束，避免DoF重排序。
"""
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
            nodes::Vector{Node}
            link::Vector{Int} # link vector of nodes       
            ptype::String # problem type: bar, beam, pstrain, pstress, sym, 3d
            dir::Int # sign of element volume, for failure checking
        end
        $name() = $name(-1,"none",Node[],Int[],"none",0)
        $name(elemtype, ptype) = $name(-1,elemtype,Node[],Int[],ptype,0)
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

function elem_x(elem::T) where T <: AbstractElem
    dim = length(elem.nodes[1].x0)
    n = length(elem.nodes)
    x = Matrix{Float64}(undef,n,dim)
    for i = 1:n, j = 1:dim
        x[i, j] = elem.nodes[i].x0[j] + elem.nodes[i].d[j]
    end
    return x
end

"""
Area of 2D elements
"""
function elem_area(elem)
    x = elem_x(elem)
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

