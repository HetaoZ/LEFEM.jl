# -----------------------------------------------
# Define an element DataType , its strain matrix,
# and integration methods.
# -----------------------------------------------

create_elem_type_and_eval(:Tet4)

"""
Tet4's strain matrix is a constant matrix.
"""
function strain_matrix(elem::Tet4, nodes, r, s)
    x = elem_x(elem, nodes)
    x1 = x[1,1]; y1 = x[1,2]
    x2 = x[2,1]; y2 = x[2,2]
    x3 = x[3,1]; y3 = x[3,2]
    x4 = x[4,1]; y4 = x[4,2]

    Nr1, Nr2, Nr3, Nr4 = r-1, 1-r, r, -r
    Ns1, Ns2, Ns3, Ns4 = s-1, -s, s, 1-s

    J11 = Nr1*x1+Nr2*x2+Nr3*x3+Nr4*x4
    J12 = Nr1*y1+Nr2*y2+Nr3*y3+Nr4*y4
    J21 = Ns1*x1+Ns2*x2+Ns3*x3+Ns4*x4
    J22 = Ns1*y1+Ns2*y2+Ns3*y3+Ns4*y4

    detJ = J11*J22-J21*J12

    b11= J22*Nr1-J12*Ns1
    b21=-J21*Nr1+J11*Ns1

    b12= J22*Nr2-J12*Ns2
    b22=-J21*Nr2+J11*Ns2

    b13= J22*Nr3-J12*Ns3
    b23=-J21*Nr3+J11*Ns3

    b14= J22*Nr4-J12*Ns4
    b24=-J21*Nr4+J11*Ns4

    B = [b11  0    b12  0    b13  0    b14  0  ;
         0    b21  0    b22  0    b23  0    b24;
         b21  b11  b22  b12  b23  b13  b24  b14] 
    return B ./ detJ, detJ
end

function elem_stress(elem::Tet4)

end

function elem_strain(elem::Tet4)

end

function elem_jacobi(elem::Tet4, nodes)
    x = elem_x(elem, nodes)
    x1 = x[1,1]; y1 = x[1,2]
    x2 = x[2,1]; y2 = x[2,2]
    x3 = x[3,1]; y3 = x[3,2]
    x4 = x[4,1]; y4 = x[4,2]

    r, s = 0.5, 0.5
    Nr1, Nr2, Nr3, Nr4 = r-1, 1-r, r, -r
    Ns1, Ns2, Ns3, Ns4 = s-1, -s, s, 1-s

    J11 = Nr1*x1+Nr2*x2+Nr3*x3+Nr4*x4
    J12 = Nr1*y1+Nr2*y2+Nr3*y3+Nr4*y4
    J21 = Ns1*x1+Ns2*x2+Ns3*x3+Ns4*x4
    J22 = Ns1*y1+Ns2*y2+Ns3*y3+Ns4*y4

    detJ = J11*J22-J21*J12
    return detJ   
end

"""
t: thickness = 1
"""
function integ_elem_elast_brick(elem::Tet4, nodes, E, ν)
    gp, gpw = GAUSS_POINT[NGP]
    D = elast_matrix_2d(E, ν, elem.ptype)
    Ke = zeros(Float64,8,8)
    for i = 1:NGP, j = 1:NGP
        B, detJ = strain_matrix(elem, nodes, gp[i], gp[j])
        Ke += (gpw[i]*gpw[j]*detJ) * (B' * D * B)
    end
    return Ke
end

"""
Accumulated mass matrix, transformed to a vector.
"""
function integ_elem_mass_brick(elem::Tet4, nodes, ρ, t)
    A = elem_area(elem, nodes)
    Me = ρ * A * t / 4 * ones(Float64,8)
    return Me
end

function integ_elem_brick(elem::Tet4, nodes, para)
    Ke = integ_elem_elast_brick(elem, nodes, para["E"], para["nu"])
    Me = integ_elem_mass_brick(elem, nodes, para["rho"], para["thickness"])
    fe = zeros(Float64, size(Me,1))
    return Ke, Me, fe
end
