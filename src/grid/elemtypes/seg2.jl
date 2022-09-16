# -----------------------------------------------
# Define an element DataType , its strain matrix,
# and integration methods.
# -----------------------------------------------

create_elem_type_and_eval(:Seg2)

"""
Seg2's strain matrix is a constant matrix.
"""
function strain_matrix(elem::Seg2, nodes)
    x = elem_x(elem,nodes)
    x1 = x[1,1]; y1 = x[1,2]
    x2 = x[2,1]; y2 = x[2,2]
    x3 = x[3,1]; y3 = x[3,2]
    D = det([1.0  x1  y1;
    1.0  x2  y2;
    1.0  x3  y3])
    b1 = y2-y3; c1 = -x2+x3
    b2 = y3-y1; c2 = -x3+x1
    b3 = y1-y2; c3 = -x1+x2
    B = [b1 0 b2 0 b3 0;
        0 c1 0 c2 0 c3;
        c1 b1 c2 b2 c3 b3]
    return B ./ D
end

function elem_stress(elem::Seg2)

end

function elem_strain(elem::Seg2)

end

function elem_jacobi(elem::Seg2, nodes)
    x = elem_x(elem, nodes)
    M = ones(Float64,3,3)
    for i = 1:3
        M[i,2:3] = x[i,:]
    end
    return det(M)    
end

"""
t: thickness = 1
"""
function integ_elem_elast_brick(elem::Seg2, nodes, E, ν)
    B = strain_matrix(elem, nodes)
    D = elast_matrix_2d(E, ν, elem.ptype)
    A = elem_area(elem, nodes)
    Ke =  B' * D * B * A
    # println("A = ",A)
    # println("B = ",B)
    # println("D = ",D)
    # println(nodes)
    return Ke
end

"""
Accumulated mass matrix, transformed to a vector.
"""
function integ_elem_mass_brick(elem::Seg2, nodes, ρ, t)
    A = elem_area(elem, nodes)
    Me = ρ * A * t / 3 * ones(Float64,6)
    return Me
end

function integ_elem_brick(elem::Seg2, nodes, para)
    Ke = integ_elem_elast_brick(elem, nodes, para["E"], para["nu"])
    Me = integ_elem_mass_brick(elem, nodes, para["rho"], para["thickness"])
    fe = zeros(Float64, size(Me,1))
    return Ke, Me, fe
end
