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

# ---------------------------
# Classes
# ---------------------------

abstract type AbstractElementType end



function create_elem_type(name)
    code = quote
        mutable struct $name <: AbstractElementType
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
# Classes
# ---------------------------
mutable struct LESystem
    K::Matrix{Float64} # stiffness
    M::Vector{Float64} # mass
    C::Matrix{Float64} # damping
    f::Vector{Float64} # force
end

# ----------------------------
# Unsupport mixed element types.
# ----------------------------
mutable struct LEStructure
    nnp::Int
    dim::Int
    ndof::Int
    nodes::Vector{Node}
    elements::Vector{T} where T <: AbstractElementType 
    boundary::Vector{Convex}
    system::LESystem
    para::Dict
    ext_f::Vector{Float64}
    cons_dof_list::Vector{Int} # constrained dofs
    cons_d_list::Vector{Real} # constrained displacements of dofs, default to be zeros.
    movable::Bool
end
LEStructure(nnp, dim, ndof, nodes, elements, boundary, system, para) = LEStructure(nnp, dim, ndof, nodes, elements, boundary, system, para, zeros(Float64, ndof), Int[], Int[], true)