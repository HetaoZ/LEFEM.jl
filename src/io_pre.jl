include("ReadGmsh.jl")
using .ReadGmsh

using DelimitedFiles
# ---------------------
# elemtype
# 1 —— Seg2
# 2 —— Tri3
# 3 —— Quad4
# ---------------------
const ET_NUMBER = Dict("Seg2"=>1, "Tri3"=>2, "Quad4"=>3)
const ET = Dict("Tri3"=>Tri3)
const ELEM_NEN = Dict(1=>2, 2=>3, 3=>4, 4=>4, 5=>8)
const ELEM_DIM = Dict("Tri3"=>2)

"""
Supported formats: .msh

elemtype: Tri3

ptype: bar, beam, pstrain, pstress, sym, 3d
"""
function read_mesh(elemtype, ptype, f::String)
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
    if nel > 0
        elements = ET[elemtype][]
        nen = ELEM_NEN[ET_NUMBER[elemtype]]
        for ie = 1:nel
            newe = ET[elemtype](elemtype, ptype)
            newe.id = ie
            newe.link = elemNodeTags[1:nen,ie]
            newe.nodes = nodes[newe.link]
            newe.dir = sign(elem_jacobi(newe))
            push!(elements, newe)
        end
        return nnp, dim, elements
    else
        error("No element detected!")
    end
end

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