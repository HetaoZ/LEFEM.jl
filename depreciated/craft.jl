"""
Integrates the polynomial in the Quad2D element. Returns the brick matrix of the element (stiffness matrix, mass matrix or damping matrix). Use default Gauss Integration Methods.
"""
function integ_elem_brick(elem::Quad2D, poly, test_poly)
    ex = poly * test_poly # integrand expression
    
end

function explicit_solver(s::LEStructure, dt)
    d = assemble_elem_field(s, :d)
    u = assemble_elem_field(s, :u)
    f = s.system.f .+ s.ext_f
    d_bk = copy(d)

    c0 = 1.0/(dt*dt)
    c1 = 0.5/dt
    c2 = 2.0*c0
    c3 = 1.0/c2

    M = diagm(s.system.M)

    a = inv(M) * (f - s.system.C*u - s.system.K*d)
    dminus = d - dt*u + c3*a

    M_eff = c0 .* M
    tmp = c0 .* M
    if s.para["damping"]
        M_eff .+= (c1 .* s.system.C)
        tmp .-= c1 .* s.system.C
    end
    f_eff = f .- (s.system.K .- c2 .* M) * d_bk .- tmp * dminus
    d = M_eff\f_eff
        
    u = (d - dminus) * c1
    a = (dminus - 2*d_bk + d) * c0
    return d, u, a
end

function newmark_beta_solver!(s::LEStructure, dt)
    d = assemble_elem_field(s, :d)
    u = assemble_elem_field(s, :u)
    a = assemble_elem_field(s, :a)
    f = s.system.f .+ s.ext_f
    u_bk = copy(u)
    a_bk = copy(a)

    δ = 0.52
    α = 0.25*(0.5+δ)^2

    K_eff = diagm(s.system.M) + α*dt^2*s.system.K
    Q = f - s.system.K * (d + dt*u + (0.5-α)*dt^2*a)
    if s.para["damping"]
        K_eff += δ*dt*s.system.C
        Q -= s.system.C * (u + (1-δ)*dt*a)
    end

    a = K_eff\Q

    u = u_bk + ((1-δ)*a_bk + δ*a) * dt
    d += u_bk*dt + ((0.5-α)*a_bk + α*a)*dt^2
    return d, u, a  
end