
# ---------------------
# elemtype
# 1 —— Seg2
# 2 —— Tri3
# 3 —— Quad4
# ---------------------
const ET_NUMBER = Dict("Seg2"=>1, "Tri3"=>2, "Quad4"=>3)
const BOUND_ET_NUMBER = Dict("Seg2"=>0, "Tri3"=>1, "Quad4"=>1) # boundary element type of a specific element type
const ET = Dict("Tri3"=>Tri3, "Quad4"=>Quad4)
const ELEM_NEN = Dict(1=>2, 2=>3, 3=>4, 4=>4, 5=>8)
const ELEM_DIM = Dict("Tri3"=>2, "Quad4"=>2)

"""
Supported formats: .msh

elemtype: Tri3

ptype: bar, beam, pstrain, pstress, sym, 3d
"""
function read_lefem_mesh(elemtype, ptype, f::String)
    dim = ELEM_DIM[elemtype]
    py_nodeTags, py_nodeCoords = get_nodes(f)
    nnp = length(py_nodeTags)
    # correct the numeration to Julia-style
    nodeTags_order = sortperm(py_nodeTags)
    nodeTags_map = Dict()
    for i = 1:nnp
        nodeTags_map[py_nodeTags[i]] = nodeTags_order[i]
    end
    nodeCoords = py_nodeCoords[1:dim,nodeTags_order]
    reordered_nodeTags = py_nodeTags[nodeTags_order] 
    nodeTags = [k for k = 1:nnp]
    # elements
    py_elemTags, py_elemNodeTags = get_elems(f, ET_NUMBER[elemtype])
    nel = length(py_elemTags)
    for i in eachindex(py_elemNodeTags)
        py_elemNodeTags[i] = nodeTags_map[py_elemNodeTags[i]]
    end
    elemNodeTags = py_elemNodeTags
    # create nodes
    nodes = tags_coords_to_nodes(nodeTags, nodeCoords)

    # elements
    if nel == 0
        error("No element detected!")
    end
    elements = ET[elemtype][]
    nen = ELEM_NEN[ET_NUMBER[elemtype]]
    for ie = 1:nel
        newe = ET[elemtype](elemtype, ptype)
        newe.id = ie
        newe.link = elemNodeTags[1:nen,ie]
        newe.dir = sign(elem_jacobi(newe, nodes))
        push!(elements, newe)
    end
    
    # boundary
    boundElemTags, boundElemNodeTags = get_bounds(f, BOUND_ET_NUMBER[elemtype], dim)


    nbound = length(boundElemTags)
    boundary = Vector{Convex}(undef, nbound)
    for i in eachindex(boundary)
        boundary[i] = Convex(boundElemTags[i], boundElemNodeTags[:,i], Float64[])
    end
    for i in eachindex(boundary)
        boundary[i].normal = outer_normal(boundary[i], nodes, get_boundary_shape!(nodes,boundary,dim), dim)
    end

    return dim, nodes, elements, boundary
end

# function relink(boundary, nodes)
#     n = length(boundary)
#     for i = 1:n
#         boundary[i] = reordernodes(boundary[i], nodes)
#     end
#     e = deepcopy(boundary)
#     for i = 2:n
#         e[i] = findnextlink(boundary, e[i-1])
#     end
#     return e
# end

# function reordernodes(face, nodes)
#     # 规定外法向量在左侧
#     if crossproduct(nodes[face.].x-b.nodes[1].x, b.n)[1] < 0 
#         tmp = deepcopy(b.nodes[1])
#         b.nodes[1] = deepcopy(b.nodes[2])
#         b.nodes[2] = tmp
#     end
#     return b
# end

# function findnextlink(edge, b)
#     for eb in edge
#         if eb.nodes[1].i == b.nodes[end].i
#             return eb
#         end
#     end
#     error("next link not found")
# end

function read_para(f::String)
    mat = readdlm(f)
    n = size(mat,1)
    para = Dict()
    for i = 1:n
        para[mat[i,1]] = mat[i,2]
    end
    return para
end

function tags_coords_to_nodes(nodeTags, nodeCoords)
    nodes = Vector{Node}(undef, length(nodeTags))
    for i = 1:length(nodeTags)
        nodes[i] = Node(size(nodeCoords,1))
        nodes[i].id = i
        nodes[i].x0 = nodeCoords[:,i]
        nodes[i].d = zeros(Float64,length(nodes[i].x0))
    end
    return nodes
end