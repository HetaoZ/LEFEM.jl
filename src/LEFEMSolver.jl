"""
Explicit
"""
module LEFEMSolver
include("LEFEMAssembly.jl")
using .LEFEMAssembly
export Tri3, LEStructure
export read_mesh, read_para, read_model, set_cons_dof!, update_system!, lefem_advance!, review
export fetch_data

using LinearAlgebra

# ---------------------------
# Common Functions
# ---------------------------

function lefem_advance!(s::LEStructure, dt, scheme)
    # read_mesh will automatically assemble the system.
    if scheme == "explicit"
        d, u, a = explicit_solver(s, dt)
    else
        error("undef")
    end
    update_elements!(s, d, u, a)
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

    # stress = s.system.K * d_bk
    # f_eff = dt^2 * (f - stress)
    # d = f_eff ./ s.system.M + 2.0*d_bk - dminus    

    
    # println("stress = ")
    # display(stress)
    # println()
    # println("f = ")
    # display(f)
    # println()
    # println("f_eff = ")
    # display(f_eff)
    # println()
    # println("d = ")
    # display(d)
    # println()
    # println("dminus = ")
    # display(dminus)
    # println()
    # println("d_bk = ")
    # display(d_bk)
    # println()
    # println("K = ")
    # display(s.system.K)
    # println()

    # println("(s.system.K .- c2 .* M)  = ")
    # display((s.system.K .- c2 .* M))
    # println()
    # println("tmp * dminus = ")
    # display(tmp * dminus)
    # println()
    # println("(s.system.K .- c2 .* M) * d_bk .- tmp * dminus = ")
    # display((s.system.K .- c2 .* M) * d_bk .- tmp * dminus)
    # println()

        
    u = (d - dminus) * c1
    a = (dminus - 2*d_bk + d) * c0
    return d, u, a
end

function update_elements!(s::LEStructure, d, u, a)
    dim = s.dim
    for ie in eachindex(s.elements)
        for k in eachindex(s.elements[ie].link)
            nid = s.elements[ie].link[k]
            s.elements[ie].nodes[k].d = d[(nid-1)*dim+1:nid*dim]
            s.elements[ie].nodes[k].u = u[(nid-1)*dim+1:nid*dim]
            s.elements[ie].nodes[k].a = a[(nid-1)*dim+1:nid*dim]
        end
    end
end


###
end
