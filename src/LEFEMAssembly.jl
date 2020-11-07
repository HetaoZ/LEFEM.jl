module LEFEMAssembly
using LinearAlgebra

include("LEFEMElem.jl")
using .LEFEMElem
export Tri3, elem_stress, elem_strain, elem_jacobi

include("io_pre.jl")
export read_mesh, read_para

# ---------------------------
# Constants
# ---------------------------
const CONSTRAIN_ALPHA = 1.0e10


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
    elements::Vector{T} where T <: AbstractElem 
    system::LESystem
    para::Dict
    ext_f::Vector{Float64}
    cons_dof_list::Vector{Int} # constrained dofs
    cons_d_list::Vector{Real} # constrained displacements of dofs, default to be zeros.
end
LEStructure(nnp, dim, ndof, elements, system, para) = LEStructure(nnp, dim, ndof, elements, system, para, zeros(Float64, ndof), Int[], Int[])
export LEStructure
# ---------------------------
# Common Functions
# ---------------------------
function assemble_K!(K, Ke, link, dim)
    dof_link = link_to_dof_link(link, dim)
    for i = 1:length(dof_link)  
        for j = 1:length(dof_link)
            K[dof_link[i], dof_link[j]] += Ke[i,j]
        end
    end 
end

function assemble_M!(M, Me, link, dim)
    dof_link = link_to_dof_link(link, dim)
    for i = 1:length(dof_link)  
        M[dof_link[i]] += Me[i]
    end     
end

function assemble_f!(f, fe, link, dim)
    dof_link = link_to_dof_link(link, dim)
    for i = 1:length(dof_link)  
        f[dof_link[i]] += fe[i]
    end 
end

function assemble!(K, M, f, Ke, Me, fe, link, dim)
    assemble_K!(K, Ke, link, dim)
    assemble_M!(M, Me, link, dim)
    assemble_f!(f, fe, link, dim)
end
export assemble!

"""
Call assemble!() only once for linear elastic problems.
"""
function assemble_system(ndof, elements, para, cons_dof_list, cons_d_list, dim)
    K = zeros(Float64, ndof, ndof)
    M = zeros(Float64, ndof)
    f = zeros(Float64, ndof)
    for elem in elements
        Ke, Me, fe = integ_elem_brick(elem, para)
        assemble!(K, M, f, Ke, Me, fe, elem.link, dim)
    end
    constrain_dof!(K, M, f, cons_dof_list, cons_d_list, dim)
    if para["damping"]
        C = diagm(para["damping coeff alpha"] .* M) .+ para["damping coeff beta"] .* K
    else
        C = zeros(Float64,1,1)
    end
    return LESystem(K, M, C, f)
end

function read_model(elemtype, ptype, meshfile, parafile)
    para = read_para(parafile)
    nnp, dim, elements = read_mesh(elemtype, ptype, meshfile)
    system = assemble_system(nnp*dim, elements, para, Int[], Int[], dim)

    s = LEStructure(nnp, dim, nnp*dim, elements, system, para)

    return s
end
export read_model

function review(s::LEStructure)
    println("-- short summary of model --")
    println("# parameters")
    for k in keys(s.para)
        println("  ",k," : ",s.para[k])
    end
    println("# mesh")
    println("  dimension : ", s.dim)
    println("  number of nodes : ", s.nnp)
    println("  number of elements : ", length(s.elements))
    println("  element type : ", s.elements[1].elemtype)
    println("  problem type : ", s.elements[1].ptype)

    println("# constrains")
    println("  number of constrained dofs : ", length(s.cons_dof_list))
    println("-- end --")
end
export review

function set_cons_dof!(s, cons_dof_list, cons_d_list)
    @assert size(cons_dof_list) == size(cons_d_list)
    for i = 1:length(cons_dof_list)
        dof = cons_dof_list[i]
        @assert dof >= 1
        @assert dof <= s.ndof
        push!(s.cons_dof_list, dof)
        push!(s.cons_d_list, cons_d_list[i])
    end
    update_system!(s)
end
export set_cons_dof!

function link_to_dof_link(link, dim)
    dof_link = Vector{Int}(undef,dim*length(link))
    for i = 1:length(link)
        for j = 1:dim
            dof_link[(i-1)*dim+j] = (link[i]-1)*dim+j
        end
    end
    return dof_link
end

function constrain_dof!(K, M, f, cons_dof_list, cons_d_list, dim)
    for k = 1:length(cons_dof_list)
        i = cons_dof_list[k]
        K[i,i] *= CONSTRAIN_ALPHA
        M[i] *= CONSTRAIN_ALPHA
        f[i] = K[i,i] * cons_d_list[k]
    end
end

function update_system!(s::LEStructure)
    s.system = assemble_system(s.ndof, s.elements, s.para, s.cons_dof_list, s.cons_d_list, s.dim)
    set_disp_of_dof!(s, s.cons_dof_list, s.cons_d_list)
end
export update_system!

function set_disp_of_dof!(s::LEStructure, cons_dof_list, cons_d_list)
    for k in eachindex(cons_dof_list)
        cons_dof = cons_dof_list[k]
        cons_d = cons_d_list[k]
        for i in eachindex(s.elements)
            dof_link = link_to_dof_link(s.elements[i].link, s.dim)
            if cons_dof in dof_link
                dof_in_elem = findall(x->x==cons_dof, dof_link)[1]
                node_k = cld(dof_in_elem, s.dim)
                disp_k = dof_in_elem%s.dim
                if disp_k == 0
                    disp_k = s.dim
                end
                s.elements[i].nodes[node_k].d[disp_k] = cons_d
            end
        end
    end
end

function assemble_elem_field(s::LEStructure, field)
    f = getfield(s.elements[1].nodes[1], field)
    d = Vector{eltype(f)}(undef, s.ndof)
    for e in s.elements
        for node in e.nodes 
            d[(node.id-1)*s.dim+1:node.id*s.dim] = getfield(node, field)
        end
    end
    return d
end
export assemble_elem_field

"""
Fetch data from elements. This differs from assemble_elem_field.
"""
function fetch_data(s::LEStructure, field)
    # f = getfield(s.elements[1].nodes[1], field)
    # d = Vector{Tuple}(undef, s.nnp)
    # for e in s.elements
    #     for node in e.nodes 
    #         d[node.id] = Tuple(getfield(node, field))
    #     end
    # end
    f = getfield(s.elements[1].nodes[1], field)
    d = Array{eltype(f)}(undef, s.dim, s.nnp)
    for e in s.elements
        for node in e.nodes 
            d[:,node.id] = getfield(node, field)
        end
    end
    return d    
    # reshape(assemble_elem_field(s,field), (s.dim, s.nnp))
end
export fetch_data
###
end