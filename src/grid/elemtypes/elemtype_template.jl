# -----------------------------------------------
# Define an element DataType , its strain matrix,
# and integration methods.
# -----------------------------------------------

create_elem_type_and_eval(:Name)

function strain_matrix(elem::Name, nodes)

end

function elem_stress(elem::Name)

end

function elem_strain(elem::Name)

end

function elem_jacobi(elem::Name, nodes)

end

function integ_elem_elast_brick(elem::Name, nodes, E, ν)
    B = strain_matrix(elem, nodes)
    D = elast_matrix_2d(E, ν, elem.ptype)
    A = elem_area(elem, nodes)
    Ke =  B' * D * B * A
    return Ke
end

"""
Accumulated mass matrix, transformed to a vector.
"""
function integ_elem_mass_brick(elem::Name, nodes, ρ, t)

end

function integ_elem_brick(elem::Name, nodes, para)
    Ke = integ_elem_elast_brick(elem, nodes, para["E"], para["nu"])
    Me = integ_elem_mass_brick(elem, nodes, para["rho"], para["thickness"])
    fe = zeros(Float64, size(Me,1))
    return Ke, Me, fe
end
