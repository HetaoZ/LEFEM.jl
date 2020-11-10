# LEFEM

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HetaoZ.github.io/LEFEM.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HetaoZ.github.io/LEFEM.jl/dev)
[![Build Status](https://github.com/HetaoZ/LEFEM.jl/workflows/CI/badge.svg)](https://github.com/HetaoZ/LEFEM.jl/actions)
[![Coverage](https://codecov.io/gh/HetaoZ/LEFEM.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HetaoZ/LEFEM.jl)

An easy-to-learn package for linear elastic dynamic FEM. This package is in active development. 

# Install
```julia
] add https://github.com/HetaoZ/LEFEM.jl.git
```

# Usage
This is an example:
```julia
# read model
s = read_model("Tri3", "pstrain", "in/rect2d.msh", "in/test_mat.para")

# constrain
set_cons_dof!(s, [1,3,5,7], [0.2,-0.2,-0.1,0.1])

# review the model info
review(s)

# set solution parameters
maxtime  = 1
maxframe = 1
cutframe = 1

t = 0
frame = 0

println(frame, "      ",t)
save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(frame))

# run the solver
while frame <= maxframe && t <= maxtime
    dt = 0.1
    lefem_advance!(s, dt, "explicit")

    global t += dt
    global frame += 1

    if frame%cutframe == 0
        println(frame, "      ",t)
        save_to_vtk(s, ["mydisp"], [:d], "out/disp_"*string(frame))
    end
end
```