
# ---------------------------
# Common Functions
# ---------------------------

function advance!(s::LEStructure, dt, scheme)
    # read_lefem_mesh will automatically assemble the system.
    if scheme == "explicit"
        d, u, a = explicit_solver(s, dt)
    else
        error("undef")
    end
    update_elements!(s, d, u, a)
    update_boundary!(s)
end

function explicit_solver(s::LEStructure, dt)
    d = assemble_elem_field(s, :d)
    u = assemble_elem_field(s, :u)
    a = assemble_elem_field(s, :a)
    f = s.system.f .+ s.ext_f
    d_bk = copy(d)

    c0 = 1.0/(dt*dt)
    c1 = 0.5/dt
    c2 = 2.0*c0
    c3 = 1.0/c2
    dminus = d - dt*u + c3*a

    M = diagm(s.system.M)

    M_eff = c0 .* M
    tmp = c0 .* M
    if s.para["damping"]
        M_eff .+= (c1 .* s.system.C)
        tmp .-= c1 .* C
    end
    f_eff = f .- (s.system.K .- c2 .* M) * d_bk .- tmp * dminus
    d = M_eff\f_eff
        
    u = (d - dminus) * c1
    a = (dminus - 2*d_bk + d) * c0
    return d, u, a
end

function update_elements!(s::LEStructure, d, u, a)
    dim = s.dim
    for nid in eachindex(s.nodes)
        s.nodes[nid].d = d[(nid-1)*dim+1:nid*dim]
        s.nodes[nid].u = u[(nid-1)*dim+1:nid*dim]
        s.nodes[nid].a = a[(nid-1)*dim+1:nid*dim]        
    end
end

function update_boundary!(s::LEStructure)
    xs = get_boundary_shape!(s)
    for k in eachindex(s.boundary)
        s.boundary[k].normal = outer_normal(s.boundary[k], s.nodes, xs, s.dim)
    end 
end

function time_step!(s::LEStructure)
    minL = 1.0
    minrho = 1.0
    for e in s.elements
        x = elem_x(e, s.nodes)
        extlink = [e.link; e.link[1]]
        for i = 1:length(e.link)
            L = norm((s.nodes[extlink[i]].x0+s.nodes[extlink[i]].d) - (s.nodes[extlink[i+1]].x0+s.nodes[extlink[i+1]].d))
            minL = min(minL, L)
        end
        minrho = min(minrho, elem_density(e, s.nodes, s.para["rho"]))
    end

    C = sqrt(s.para["E"]/minrho)
    dt = pi * minL / C
    return dt
end