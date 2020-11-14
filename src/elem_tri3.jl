# -----------------------------------------------
# Define an element DataType , its strain matrix,
# and integration methods.
# -----------------------------------------------

create_elem_type_and_eval(:Tri3)

"""
Tri3's strain matrix is a constant matrix.
"""
function strain_matrix(elem::Tri3, nodes)
    x = elem_x(elem,nodes)
    M = ones(Float64,3,3)
    for i = 1:3
        M[i,2:3] = x[i,:]
    end
    D = det(M)
    # B = zeros(Float64,3,6)
    # for k = 1:3
    #     p = Bool[true,true,true]
    #     p[k] = false
    #     MP = M[p,:]
    #     b = -det(MP[:,[1,3]])
    #     c = det(MP[:,[1,2]])
    #     Bk = zeros(Float64,3,2)
    #     Bk[1,1] = b
    #     Bk[2,2] = c
    #     Bk[3,1] = c
    #     Bk[3,2] = b
    #     B[:, 2*k-1:2*k] = Bk
    # end
    x1 = x[1,1]; y1 = x[1,2]
    x2 = x[2,1]; y2 = x[2,2]
    x3 = x[3,1]; y3 = x[3,2]
    b1 = y2-y3; c1 = -x2+x3
    b2 = y3-y1; c2 = -x3+x1
    b3 = y1-y2; c3 = -x1+x2
    B = [b1 0 b2 0 b3 0;
        0 c1 0 c2 0 c3;
        c1 b1 c2 b2 c3 b3]
    return B ./ D
end

function elem_stress(elem::Tri3)

end

function elem_strain(elem::Tri3)

end

function elem_jacobi(elem::Tri3, nodes)
    x = elem_x(elem, nodes)
    M = ones(Float64,3,3)
    for i = 1:3
        M[i,2:3] = x[i,:]
    end
    return det(M)    
end

"""
t: thickness
"""
function integ_elem_elast_brick(elem::Tri3, nodes, E, ν)
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
function integ_elem_mass_brick(elem::Tri3, nodes, ρ, t)
    A = elem_area(elem, nodes)
    Me = ρ * A * t / 3 * ones(Float64,6)
    return Me
end

function integ_elem_brick(elem::Tri3, nodes, para)
    Ke = integ_elem_elast_brick(elem, nodes, para["E"], para["nu"])
    Me = integ_elem_mass_brick(elem, nodes, para["rho"], para["thickness"])
    fe = zeros(Float64, size(Me,1))
    return Ke, Me, fe
end
