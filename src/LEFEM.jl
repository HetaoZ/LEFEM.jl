module LEFEM
using LinearAlgebra, DelimitedFiles, WriteVTK, ReadGmsh, MathKits, Statistics, PointInPoly

include("base.jl")
include("elem.jl")
include("io_pre.jl")
include("assembly.jl")
export read_model, review, set_cons_dof!, fetch_data, LEStructure
include("solver.jl")
export lefem_advance!
include("io_post.jl")
export save_to_vtk

###
end