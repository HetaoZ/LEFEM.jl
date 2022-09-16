
# ---------------------------
# Constants
# ---------------------------


# ---------------------------
# Material definitions
# ---------------------------
include("grid/material/linear_elasticity.jl")

# ---------------------------
# Element definitions
# ---------------------------
include("grid/elemtypes/tri3.jl")
include("grid/elemtypes/quad4.jl")
# include("elemtypes/tet4.jl")
# include("elemtypes/hex8.jl")
# ---------------------------
# Common Functions
# ---------------------------

function elem_x(elem::T, nodes::Array{Node}) where T <: AbstractElementType
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


function polygon_area(x, y)
    n = length(x)
    if n < 3 
        return 0.
    end
    s = x[n]*(y[1] - y[n-1]) + x[1]*(y[2] - y[n])
    for i = 2:n-1
        s += x[i]*(y[i+1] - y[i-1])
    end
    return s * 0.5
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



