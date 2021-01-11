module LEFEM
using LinearAlgebra, DelimitedFiles, WriteVTK, ReadGmsh, MathKits, Statistics, PointInPoly

const GAUSS_POINT = Dict(
                         1 => (0.5, 1.0),
                         2 => ([0.21132486540518702,  0.788675134594813], [0.5, 0.5]),
                         3 => ([0.11270166537925852, 0.5, 0.8872983346207415],  [0.277777777777778, 0.4444444444444445, 0.277777777777778])
)  
const NGP = 2

include("base.jl")
include("elem.jl")
include("copy.jl")
include("io_pre.jl")
include("assembly.jl")
export read_model, review, cons_dof!, cons_dof_in_box!, cons_force_in_box!, fetch_data, LEStructure
include("solver.jl")
export advance!, time_step!
include("io_post.jl")
export save_to_vtk

###
end