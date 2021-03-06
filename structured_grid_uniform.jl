include("./module.jl")
include("./controls.jl")
using .mesh

function structured_grid_uniform!(
    ๐::controls,
    cell::Vector{mesh.Cell},
    face::Vector{mesh.Face},
    face_internal::Vector{mesh.Face},
    face_boundary::Vector{mesh.Face},
    face_boundary_top::Vector{mesh.Face},
    face_boundary_bottom::Vector{mesh.Face},
    face_boundary_left::Vector{mesh.Face},
    face_boundary_right::Vector{mesh.Face}
    )

    N = ๐.Nx*๐.Ny*๐.Nz

    ฮx = ๐.Lx/๐.Nx
    ฮy = ๐.Ly/๐.Ny 
    ฮz = ๐.Lz/๐.Nz 

    x = 0.5*ฮx:ฮx:๐.Lx
    y = 0.5*ฮy:ฮy:๐.Ly
    z = 0.5*ฮz:ฮz:๐.Lz

    length(y)

    for i in 1:N 
        push!(cell, mesh.Cell())
    end

    #println(N)
    #println(length(cell))


    for i in 1:๐.Nx 
        for j in 1:๐.Ny
            k=1
            ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
            cell[ijk].x = x[i]
            cell[ijk].y = y[j]
            cell[ijk].z = 0.5*ฮz #z[i]

            cell[ijk].ฮฉ = ฮx*ฮy*ฮz

            cell[ijk].Qแต = zeros(Float64,6)
            cell[ijk].Qโฟ = zeros(Float64,6)
            cell[ijk].Qโฟโปยน = zeros(Float64,6)
        end
    end



    #println(length(cell))
    # internal faces
    for i in 2:๐.Nx
        for j in 2:๐.Ny
            k = 1
            ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
            imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
            ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
            nฬ = [0.0, 1.0, 0.0]
            ฮS = ฮx*ฮz
            owner = ijmk
            neighbour = ijk
            push!(face, mesh.Face(cell[ijk].x, ฮy*(j-1), 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
            nฬ = [1.0, 0.0, 0.0]
            ฮS = ฮy*ฮz
            owner = imjk
            neighbour = ijk
            push!(face, mesh.Face(ฮx*(i-1), cell[ijk].y, 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
        end
    end

    for j in 2:๐.Ny
        i = 1
        k = 1
        ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
        nฬ = [0.0, 1.0, 0.0]
        ฮS = ฮx*ฮz
        owner = ijmk
        neighbour = ijk
        push!(face, mesh.Face(cell[ijk].x, ฮy*(j-1), 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
    end

    for i in 2:๐.Nx
        j = 1
        k = 1
        ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
        nฬ = [1.0, 0.0, 0.0]
        ฮS = ฮy*ฮz
        owner = imjk
        neighbour = ijk
        push!(face, mesh.Face(ฮx*(i-1), cell[ijk].y, 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
    end

    face_internal_num = length(face)
    for i in face
        push!(face_internal, i)
    end

    # boundary faces
    # Left
    for j in 1:๐.Ny
        i = 1
        k = 1
        ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
        nฬ = [-1.0, 0.0, 0.0]
        ฮS = ฮy*ฮz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(ฮx*(i-1), cell[ijk].y, 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
    end

    face_total_left = length(face)
    for i in face_internal_num+1:face_total_left
        push!(face_boundary_left, face[i])
    end

    # Bottom
    for i in 1:๐.Nx
        j = 1
        k = 1
        ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
        nฬ = [0.0, -1.0, 0.0]
        ฮS = ฮx*ฮz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(cell[ijk].x, ฮy*(j-1), 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
    end
    
    face_total_bottom = length(face)
    for i in face_total_left+1:face_total_bottom
        push!(face_boundary_bottom, face[i])
    end

    # Right
    for j in 1:๐.Ny
        i = ๐.Nx
        k = 1
        ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
        nฬ = [1.0, 0.0, 0.0]
        ฮS = ฮy*ฮz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(ฮx*๐.Nx, cell[ijk].y, 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
    end
    
    face_total_right = length(face)
    for i in face_total_bottom+1:face_total_right
        push!(face_boundary_right, face[i])
    end

    # Top
    for i in 1:๐.Nx
        j = ๐.Ny
        k = 1
        ijk = i + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        imjk = (i-1) + ๐.Nx*(j-1) + ๐.Nx*๐.Ny*(k-1)
        ijmk = i + ๐.Nx*(j-2) + ๐.Nx*๐.Ny*(k-1)
        nฬ = [0.0, 1.0, 0.0]
        ฮS = ฮx*ฮz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(cell[ijk].x, ฮy*๐.Ny, 0.5*ฮz, owner, neighbour, nฬ, ฮS, [], [], []))
    end
    
    face_total_top = length(face)
    for i in face_total_right+1:face_total_top
        push!(face_boundary_top, face[i])
    end

    
    face_total = length(face)
    for i in face_internal_num+1:face_total
        push!(face_boundary, face[i])
    end
    

    for i in cell
        for j in 1:30
            push!(i.var,0.0)
        end
    end

    for i in face
        for j in 1:9
            push!(i.varโ,0.0)
            push!(i.varแตฃ,0.0)
        end
    end
    #=
    =#

    return nothing

end