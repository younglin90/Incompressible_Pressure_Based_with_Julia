

function momentum!(
    👉::controls,
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face},
    faces_boundary_top::Vector{mesh.Face},
    faces_boundary_bottom::Vector{mesh.Face},
    faces_boundary_left::Vector{mesh.Face},
    faces_boundary_right::Vector{mesh.Face}
    )

    A_rows::Vector{Int64} = []
    A_cols::Vector{Int64} = []
    A_vals::Vector{Float64} = []
    
    # contruct A matrix diagonal terms
    # contruct B vector
    B = zeros(Float64, length(cells), 2)
    
    diagon = 1

    for cell in cells
        
        push!(A_rows, diagon)
        push!(A_cols, diagon)

        tmp_A_var = cell.var[👉.ρ]*cell.Ω/👉.Δt

        push!(A_vals, tmp_A_var)

        tmp_Bx_var = -( 
            cell.var[👉.ρ]*cell.var[👉.u] - cell.var[👉.ρⁿ]*cell.var[👉.uⁿ]
            )*cell.Ω/👉.Δt
            
        tmp_By_var = -( 
            cell.var[👉.ρ]*cell.var[👉.v] - cell.var[👉.ρⁿ]*cell.var[👉.vⁿ]
            )*cell.Ω/👉.Δt

        tmp_By_var += cell.var[👉.ρ]*cell.Ω * (-9.8)
        #tmp_By_var += cell.var[👉.ρⁿ]*cell.Ω * (-9.8)

        B[diagon, 1] = tmp_Bx_var
        B[diagon, 2] = tmp_By_var

        diagon += 1

    end

    
    ∂Δp∂x0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pₙ = 0.5 * (cells[face.owner].var[👉.p] + cells[face.neighbour].var[👉.p])
        ∂Δp∂x0[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.neighbour, 1] -= pₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂Δp∂x0[face.neighbour, 2] -= pₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂Δp∂x0[face.neighbour, 3] -= pₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end

    for face in faces_boundary
        pₙ = cells[face.owner].var[👉.p]
        ∂Δp∂x0[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end

    # contruct A matrix  
    # contruct B vector 
    for face in faces_internal

        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        μₗ = cells[face.owner].var[👉.μ]
        μᵣ = cells[face.neighbour].var[👉.μ]
        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        #invρΔt = (wₗ/ρₗ + wᵣ/ρᵣ) * 👉.Δt
        invρΔt = 0.5 * (1.0/ρₗ + 1.0/ρᵣ) * 👉.Δt
        
        # Rhie-Chow
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += 0.5 * 👉.Δt / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        #=
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += invρΔt * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += invρΔt * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        =#
        Uₙ -= invρΔt * (pᵣ-pₗ) / ΔLR

        wₗ = sign(Uₙ)
        wᵣ = 1.0 - wₗ

        ρₙ = wₗ * ρₗ + wᵣ * ρᵣ
        uₙ = wₗ * uₗ + wᵣ * uᵣ
        vₙ = wₗ * vₗ + wᵣ * vᵣ
        μₙ = wₗ * μₗ + wᵣ * μᵣ

        push!(A_rows, face.owner)
        push!(A_cols, face.neighbour)

        tmp_AL_var = wᵣ * ρᵣ * Uₙ * ΔS
        #tmp_AL_var -= μₙ / ΔLR * ΔS
        push!(A_vals, tmp_AL_var)
        
        push!(A_rows, face.neighbour)
        push!(A_cols, face.owner)

        tmp_AR_var = -wₗ * ρₗ * Uₙ * ΔS 
        #tmp_AR_var -= μₙ / ΔLR * ΔS
        push!(A_vals, tmp_AR_var)

        A_vals[face.owner] += wₗ * ρₗ * Uₙ * ΔS 
        A_vals[face.neighbour] -= wᵣ * ρᵣ * Uₙ * ΔS 

        # convective terms
        #B[face.owner, 1] -= ρₗ * uₙ * Uₙ * ΔS
        #B[face.neighbour, 1] += ρᵣ * uₙ * Uₙ * ΔS

        #B[face.owner, 2] -= ρₗ * vₙ * Uₙ * ΔS
        #B[face.neighbour, 2] += ρᵣ * vₙ * Uₙ * ΔS

        covflux = ( wₗ * ρₗ * uₗ + wᵣ * ρᵣ * uᵣ ) * Uₙ * ΔS
        B[face.owner, 1] -= covflux
        B[face.neighbour, 1] += covflux

        covflux = ( wₗ * ρₗ * vₗ + wᵣ * ρᵣ * vᵣ ) * Uₙ * ΔS
        B[face.owner, 2] -= covflux
        B[face.neighbour, 2] += covflux

        # pressure terms
        pₙ = 0.5 * (pₗ + pᵣ)

        B[face.owner, 1] -= pₙ * face.n̂[1] * ΔS
        B[face.neighbour, 1] += pₙ * face.n̂[1] * ΔS 
        
        B[face.owner, 2] -= pₙ * face.n̂[2] * ΔS
        B[face.neighbour, 2] += pₙ * face.n̂[2] * ΔS

#=
        # viscous terms
        B[face.owner, 1] += μₙ * (uᵣ - uₗ) / ΔLR * ΔS
        B[face.neighbour, 1] -= μₙ * (uᵣ - uₗ) / ΔLR * ΔS
        
        B[face.owner, 2] += μₙ * (vᵣ - vₗ) / ΔLR * ΔS
        B[face.neighbour, 2] -= μₙ * (vᵣ - vₗ) / ΔLR * ΔS
=#

    end
    


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ρₗ = cells[face.owner].var[👉.ρ]
        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]
        Uₙ0 = Uₙ

        invU = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = invU * face.n̂[1]
        Uₙ += invV * face.n̂[2]
        Uₙ += invW * face.n̂[3]
        
        #A_vals[face.owner] += ρₗ * Uₙ * ΔS

        Uₙ = 0.0

        # convective terms
        B[face.owner, 1] -= ρₗ * invU * Uₙ * ΔS
        B[face.owner, 2] -= ρₗ * invV * Uₙ * ΔS
        
        # pressure terms
        pₙ = cells[face.owner].var[👉.p]
        B[face.owner, 1] -= pₙ * face.n̂[1] * ΔS
        B[face.owner, 2] -= pₙ * face.n̂[2] * ΔS
        
    end
 
    A = sparse(A_rows,A_cols,A_vals)
    #x = zeros(Float64, length(cells), 1)

    #println(A)
    #println(B)

    #spy(A, marker=".", markersize=1)
    #gui()

    ps = MKLPardisoSolver()
    ΔU = solve(ps, A, B)
    
#    P = ilu(A, τ = 0.1)

#    Δu = gmres!(x, A, Bx, Pl = P, log=true, maxiter = 1000)
#    Δv = gmres!(x, A, By, Pl = P, log=true, maxiter = 1000)
    #println(maximum(Δu))
    
    diagon = 1
    for cell in cells

        cell.var[👉.u] += 0.7*ΔU[diagon, 1]
        cell.var[👉.v] += 0.7*ΔU[diagon, 2]

        diagon += 1
    end


    return log10(norm(ΔU))
   

end
