using WriteVTK

const VTK_CELL_TYPE = Dict("Tri3"=>VTKCellTypes.VTK_TRIANGLE)

function create_vtkfile(s::LEStructure, fname)
    points, cells = pts_cells(s)
    f = vtk_grid(fname, points, cells)
    return f
end

function pts_cells(s::LEStructure)
    pts = Matrix{Float64}(undef, s.dim, s.nnp)
    cells = Vector{MeshCell}(undef, length(s.elements))
    for eid = 1:length(s.elements)
        e = s.elements[eid]
        for k = 1:length(e.link)
            node = s.nodes[e.link[k]]
            pts[:, node.id] = node.x0
        end
        cells[eid] = MeshCell(VTK_CELL_TYPE[e.elemtype], e.link)
    end
    return pts, cells
end

function save_to_vtk(s, datanames, fields, fname)
    vtkfile = create_vtkfile(s, fname)
    for i in eachindex(datanames)
        vtkfile[datanames[i]] = fetch_data(s, fields[i])
    end
    outfiles = vtk_save(vtkfile)
    println("       saved to ",outfiles[1])
end
