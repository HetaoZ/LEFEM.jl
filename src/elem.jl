
# ---------------------------
# Constants
# ---------------------------


# ---------------------------
# Material definitions
# ---------------------------
include("material.jl")

# ---------------------------
# Element definitions
# ---------------------------
include("elemtypes/tri3.jl")
include("elemtypes/quad4.jl")
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
function elem_area(elem, nodes)
    x = elem_x(elem,nodes)
    A = polygon_area(x[:,1], x[:,2])
    return A
end

function elem_density(elem, nodes, rho0)
    dim = length(nodes[1].x0)
    n = length(elem.link)
    x = Matrix{Float64}(undef,n,dim)
    for i = 1:n, j = 1:dim
        x[i, j] = nodes[elem.link[i]].x0[j]
    end
    A = polygon_area(x[:,1], x[:,2])
    mass = A * rho0
    rho = mass / elem_area(elem, nodes)
    return rho
end

"""
Check if an element fails.
"""
function check_elem(elem)
    if sign(elem_jacobi(elem)) != sign(elem.dir)
        error("Element ( "*string(elem.id)*" ) failed.")
    end
end



