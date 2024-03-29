
# ---------------------------
# Constants
# ---------------------------
const CONSTRAIN_ALPHA = 1.0e10

# ---------------------------
# Common Functions
# ---------------------------
function assemble_K!(K, Ke, link, dim)
    dof_link = link_to_dof_link(link, dim)
    for i = eachindex(dof_link)  
        for j = eachindex(dof_link)
            K[dof_link[i], dof_link[j]] += Ke[i,j]
        end
    end 
end

function assemble_M!(M, Me, link, dim)
    dof_link = link_to_dof_link(link, dim)
    for i = eachindex(dof_link)  
        M[dof_link[i]] += Me[i]
    end     
end

function assemble_f!(f, fe, link, dim)
    dof_link = link_to_dof_link(link, dim)
    for i = eachindex(dof_link)  
        f[dof_link[i]] += fe[i]
    end 
end

function assemble!(K, M, f, Ke, Me, fe, link, dim)
    assemble_K!(K, Ke, link, dim)
    assemble_M!(M, Me, link, dim)
    assemble_f!(f, fe, link, dim)
end

"""
Call assemble!() only once for linear elastic problems.
"""
function assemble_system(ndof, nodes, elements, para, cons_dof_list, cons_d_list, dim)
    K = zeros(Float64, ndof, ndof)
    M = zeros(Float64, ndof)
    f = zeros(Float64, ndof)
    for elem in elements
        Ke, Me, fe = integ_elem_brick(elem, nodes, para)
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

"""
Sometimes it returns `ERROR: gmshModelGetBoundary returned non-zero error code: 1`. But there should be no error. Just call it again.
"""
function read_model(elemtype, ptype, meshfile, parafile)
    try
        para = read_para(parafile)
        dim, nodes, elements, boundary = read_lefem_mesh(elemtype, ptype, meshfile)
        nnp = length(nodes)
        system = assemble_system(nnp*dim, nodes, elements, para, Int[], Int[], dim)
    
        s = LEStructure(nnp, dim, nnp*dim, nodes, elements, boundary, system, para)
    
        return s
    catch
        para = read_para(parafile)
        dim, nodes, elements, boundary = read_lefem_mesh(elemtype, ptype, meshfile)
        nnp = length(nodes)
        system = assemble_system(nnp*dim, nodes, elements, para, Int[], Int[], dim)
    
        s = LEStructure(nnp, dim, nnp*dim, nodes, elements, boundary, system, para)
    
        return s
    end        
end

function review(s::LEStructure)
    println("-- short summary of LEStructure --")
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
    println()
end

function cons_dof!(s, cons_dof_list, cons_d_list)
    @assert size(cons_dof_list) == size(cons_d_list)
    for i = eachindex(cons_dof_list)
        dof = cons_dof_list[i]
        @assert dof >= 1
        @assert dof <= s.ndof
        push!(s.cons_dof_list, dof)
        push!(s.cons_d_list, cons_d_list[i])
    end
    update_system!(s)
end

@inline function betweeneq(a::Vector{T}, lo::Vector{T}, hi::Vector{T}) where T <: Real
    return all(lo .≤ a .≤ hi)
end

function cons_dof_in_box!(s, point1, point2; axis = "all", d = "zero")
    for node in s.nodes
        if betweeneq(node.x0 + node.d, point1, point2)
            if axis == "all"
                # Constrain dofs in all axis.
                dof = [s.dim*(node.id-1)+k for k=1:s.dim]
                if d == "zero"
                    cons_d = [0 for k=1:s.dim]
                else
                    @assert length(d) == s.dim
                    cons_d = d
                end
            elseif axis in (1,2,3) 
                dof = [s.dim*(node.id-1)+axis]
                if d == "zero"
                    cons_d = [0]
                else
                    @assert typeof(d) <: Real
                    cons_d = [d]
                end
            else
                error("undefined axis")
            end
            
            append!(s.cons_dof_list, dof)
            append!(s.cons_d_list, cons_d)
        end
    end
    update_system!(s)
end

function cons_force_in_box!(s, point1, point2, force)
    @assert length(force) == s.dim
    for node in s.nodes
        if betweeneq(node.x0 + node.d, point1, point2)
            s.ext_f[s.dim*(node.id-1)+1:s.dim*node.id] = force
        end
    end
end

function link_to_dof_link(link, dim)
    dof_link = Vector{Int}(undef,dim*length(link))
    for i = eachindex(link)
        for j = 1:dim
            dof_link[(i-1)*dim+j] = (link[i]-1)*dim+j
        end
    end
    return dof_link
end

function constrain_dof!(K, M, f, cons_dof_list, cons_d_list, dim)
    for k = eachindex(cons_dof_list)
        i = cons_dof_list[k]
        K[i,i] *= CONSTRAIN_ALPHA
        M[i] *= CONSTRAIN_ALPHA
        f[i] = K[i,i] * cons_d_list[k]
    end
end

function update_system!(s::LEStructure)
    s.system = assemble_system(s.ndof, s.nodes,s.elements, s.para, s.cons_dof_list, s.cons_d_list, s.dim)
    set_disp_of_dof!(s, s.cons_dof_list, s.cons_d_list)
end

function set_disp_of_dof!(s::LEStructure, cons_dof_list, cons_d_list)
    for k in eachindex(cons_dof_list)
        cons_dof = cons_dof_list[k]
        cons_d = cons_d_list[k]
        cons_node = ceil(Int, cons_dof / s.dim)
        cons_node_k = cons_dof - (cons_node-1)*s.dim
        s.nodes[cons_node].d[cons_node_k] = cons_d
    end
end

function assemble_elem_field(s::LEStructure, field)
    f = getfield(s.nodes[1], field)
    d = Vector{eltype(f)}(undef, s.ndof)
    for node in s.nodes 
        d[(node.id-1)*s.dim+1:node.id*s.dim] = getfield(node, field)
    end
    return d
end

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
    f = getfield(s.nodes[1], field)
    d = Array{eltype(f)}(undef, s.dim, s.nnp)
    for node in s.nodes 
        d[:,node.id] = getfield(node, field)
    end
    return d    
    # reshape(assemble_elem_field(s,field), (s.dim, s.nnp))
end

# --------------------
# assemble the `edge`
# --------------------

# function relink(edge)
#     n = length(edge)
#     for i = 1:n
#         edge[i] = reordernodes(edge[i])
#     end
#     e = deepcopy(edge)
#     for i = 2:n
#         e[i] = findnextlink(edge, e[i-1])
#     end
#     return e
# end

# function reordernodes(b)
#     # 规定法向量在左侧
#     if crossproduct(b.nodes[2].x-b.nodes[1].x, b.n)[1] < 0 
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

function get_boundary_shape!(s::LEStructure)
    if s.dim == 2
        n = length(s.boundary)
        x = Matrix{Float64}(undef, n, s.dim)
        for k = 1:n
            c = s.boundary[k]
            node = s.nodes[c.link[1]]
            x[k, :] = node.x0 + node.d
        end
        return x
    else
        error("undef")
    end
end

function get_boundary_shape!(nodes, boundary, dim)
    if dim == 2
        n = length(boundary)
        x = Matrix{Float64}(undef, n, dim)
        for k = 1:n
            c = boundary[k]
            node = nodes[c.link[1]]
            x[k, :] = node.x0 + node.d
        end
        return x
    else
        error("undef")
    end
end

function link_to_faces(c::Convex)
    return ntuple(i->(c.link[i], c.link[i+1]), length(c.link)-1)
end

function outer_normal(c::Convex, nodes::Vector{Node}, xs::Array{Float64}, dim::Int)
    normal = convex_normal(c, nodes, dim)
    nodesx = map(k->nodes[k].x0+nodes[k].d, c.link)
    xc = sum(nodesx)/length(nodesx)
    bias = 1.e-10
    xt = xc + normal*bias
    if PointInPoly.pinpoly(ntuple(i->Tuple(xs[:,i]),size(xs,2))..., link_to_faces(c), Tuple(xt)) == 1
        normal *= -1
    end
    return normal
end

function rotate_matrix(angle::Real)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end

function convex_normal(c, nodes, dim)
    if dim == 1
        normal = [1]
    elseif dim == 2
        normal = rotate_matrix(pi/2) * ((nodes[c.link[2]].x0+nodes[c.link[2]].d) - ((nodes[c.link[1]].x0+nodes[c.link[1]].d)))
    else
        error("undef")
    end
    return truncated_normalize(normal)
end

function truncated_normalize(v)
    v = normalize(v)
    a = map(x->x^2, v)
    p = sortperm(v)
    for i = eachindex(v)
        v[p[i]] = sqrt(sum(a) - sum(a[p[[(k==i ? false : true) for k=eachindex(v)]]]))
    end
    return v
end