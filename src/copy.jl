function copy_structure!(s::LEStructure)
    return LEStructure(s.nnp,
            s.dim,
            s.ndof,
            copy_nodes!(s.nodes),
            copy_elements!(s.elements),
            copy_boundary!(s.boundary),
            deepcopy(s.system),
            s.para,
            copy(s.ext_f),
            copy(s.cons_dof_list),
            copy(s.cons_d_list),
            s.movable)
end

function copy_nodes!(nodes::Vector{Node})
    return deepcopy.(nodes)
end

function copy_elements!(elements::Vector{T} where T <: AbstractElementType )
    return deepcopy.(elements)
end

function copy_boundary!(boundary::Vector{Convex})
    return deepcopy.(boundary)
end