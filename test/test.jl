try 
    using .LEFEM
catch
    include("../src/LEFEM.jl")
    using .LEFEM
end

# read model
s = read_model("Tri3", "pstrain", "in/plate.msh", "in/steel.para")

# constrain
cons_dof_in_box!(s, [-1,-1e-7], [2,1e-7])
cons_force_in_box!(s, [-1e-7,0.02-1e-7], [1e-7,0.02+1e-7], [1e3, 0])

# # review the model info
review(s)

# set solution parameters
maxtime  = 1
maxframe = 10
cutframe = 1
N = 1000000

t = 0
frame = 0

println(frame, "      ",t)
save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(N+frame))

# run the solver
while frame <= maxframe && t <= maxtime
    dt = 1e-6
    advance!(s, dt, "explicit")

    global t += dt
    global frame += 1

    # if frame%cutframe == 0
        println(frame, "      ",t)
        save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(N+frame))

        display(fetch_data(s,:d)); println()
        display(fetch_data(s,:u)); println()
        display(fetch_data(s,:a)); println()
        println("-----------------")
    # end
end