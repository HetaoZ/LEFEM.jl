include("../src/LEFEM.jl")
using .LEFEM

# read model
s = read_model("Quad4", "pstrain", "in/simple2d.msh", "in/test_mat.para")

# constrain
cons_dof_in_box!(s, [0.5,-10], [1.5,10])
cons_force_in_box!(s, [-1e-7,-1e-7], [1e-7,1e-7], [1e6, 0])

# # review the model info
review(s)

# set solution parameters
maxtime  = 1
maxframe = 1
N = 1000000

t = 0
frame = 0

println(frame, "      ",t)
save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(N+frame))

# run the solver
while frame <= maxframe && t <= maxtime
    dt = 1e-6
    advance!(s, dt, "newmark")

    global t += dt
    global frame += 1

    if frame%20 == 0
        println(frame, "      ",t)
        save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(N+frame))
    end
end