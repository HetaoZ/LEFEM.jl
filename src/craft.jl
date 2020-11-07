"""
Integrates the polynomial in the Quad2D element. Returns the brick matrix of the element (stiffness matrix, mass matrix or damping matrix). Use default Gauss Integration Methods.
"""
function integ_elem_brick(elem::Quad2D, poly, test_poly)
    ex = poly * test_poly # integrand expression
    
end