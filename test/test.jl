try 
    using .LEFEM
catch
    include("../src/LEFEM.jl")
    using .LEFEM
end

# read model
s = read_model("Tri3", "pstrain", "in/rect2d.msh", "in/test_mat.para")

# constrain
set_cons_dof!(s, [1,3,5,7], [0.2,-0.2,-0.1,0.1])

# "乘大数法。但该方法用于动力学时还需考虑M和C矩阵的变化，所以这是有bug的。对M和C的元素也乘以大数即可解决。"

# review the model info
review(s)

# set solution parameters
maxtime  = 1
maxframe = 10
cutframe = 1

t = 0
frame = 0

print(frame, "      ",t)
save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(frame))
println()

# run the solver
while frame <= maxframe && t <= maxtime
    dt = 0.1
    lefem_advance!(s, dt, "explicit")

    global t += dt
    global frame += 1

    if frame%cutframe == 0
        print(frame, "      ",t)
        save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(frame))
    end
end