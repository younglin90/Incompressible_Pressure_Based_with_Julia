

function coupled!(
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

    B_n = 3
    A_n = B_n * B_n

    A_rows = zeros(Int64, length(cells)*A_n)
    A_cols = zeros(Int64, length(cells)*A_n)
    A_vals = zeros(Float64, length(cells)*A_n)
    B = zeros(Float64, length(cells)*B_n)
    
    diagon = 1

    for cell in cells
        
        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        Ω = cell.Ω
        Δt = 👉.Δt
        u = cell.var[👉.u]
        v = cell.var[👉.v]
        ρ = cell.var[👉.ρ]
        p = cell.var[👉.p]
        Y₁ = cell.var[👉.Y₁]
        ρⁿ = cell.var[👉.ρⁿ]
        uⁿ = cell.var[👉.uⁿ]
        vⁿ = cell.var[👉.vⁿ]
        pⁿ = cell.var[👉.pⁿ]
        Y₁ⁿ = cell.var[👉.Y₁ⁿ]

        # continuity
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 1
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 1; A_cols[i] = ijStart + 3
        
        B[ijStart + 1] = 0.0

        # x-momentum
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 1

        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 2
        A_vals[i] = ρ*Ω/Δt
        
        i += 1
        A_rows[i] = ijStart + 2; A_cols[i] = ijStart + 3


        B[ijStart + 2] = -ρ*(u - uⁿ)*Ω / Δt

        # y-momentum
        g = -9.8
        #g = 0.0

        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 1
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 2
        
        i += 1
        A_rows[i] = ijStart + 3; A_cols[i] = ijStart + 3
        A_vals[i] = ρ*Ω/Δt


        B[ijStart + 3] = -ρ*(v - vⁿ)*Ω / Δt + ρ*g*Ω 



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
        
        ijStartₗ = B_n*(face.owner-1)
        ijStartᵣ = B_n*(face.neighbour-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        wₗ = 0.0#cells[face.owner].var[👉.w]
        wᵣ = 0.0#cells[face.neighbour].var[👉.w]
        Y₁ₗ = cells[face.owner].var[👉.Y₁]
        Y₁ᵣ = cells[face.neighbour].var[👉.Y₁]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        ΔS = face.ΔS

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        d̂ = 👉.Δt / ρˢ
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR

        
        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        wₙ = 0.0#Wₗ * wₗ + Wᵣ * wᵣ
        Y₁ₙ = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ

        pₙ = 0.5 * (pₗ + pᵣ)

        
        iₗ = A_n*(face.owner-1)
        iᵣ = A_n*(face.neighbour-1)



        #------------------------
        # continuity
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( - d̂ / ΔLR * ΔS ))
        
        A_vals[iᵣ] -= ( - d̂ / ΔLR * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( d̂ / ΔLR * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        
        A_vals[iₗ] += ( 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( 0.5 * face.n̂[1] * ΔS ))
        
        A_vals[iᵣ] -= ( 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( 0.5 * face.n̂[2] * ΔS ))
        
        A_vals[iᵣ] -= ( 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( 0.5 * face.n̂[2] * ΔS ))
        
        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( 0.5 * face.n̂[1] * ΔS + ρₗ * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( 0.5 * face.n̂[1] * ΔS - ρₗ * d̂ / ΔLR * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( 0.5 * face.n̂[1] * ΔS - ρᵣ * d̂ / ΔLR * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( 0.5 * face.n̂[1] * ΔS + ρᵣ * d̂ / ΔLR * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₗ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₗ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρₗ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρᵣ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρᵣ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρᵣ * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1
        A_vals[iₗ] += ( ρₗ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₗ * 0.5 * face.n̂[2] * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρᵣ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρᵣ * 0.5 * face.n̂[2] * uₙ * ΔS ))


        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] +=  ( 0.5 * face.n̂[2] * ΔS + ρₗ * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( 0.5 * face.n̂[2] * ΔS - ρₗ * d̂ / ΔLR * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( 0.5 * face.n̂[2] * ΔS - ρᵣ * d̂ / ΔLR * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( 0.5 * face.n̂[2] * ΔS + ρᵣ * d̂ / ΔLR * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₗ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₗ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρᵣ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρᵣ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₗ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₗ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρₗ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρᵣ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρᵣ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρᵣ * Uₙ * ΔS ))

        # ----------------------------

        # B
        B[ijStartₗ + 1] -= ( Uₙ * ΔS )
        B[ijStartᵣ + 1] += ( Uₙ * ΔS )
        
        B[ijStartₗ + 2] -= ( ρₗ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartᵣ + 2] += ( ρᵣ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )

        B[ijStartₗ + 3] -= ( ρₗ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartᵣ + 3] += ( ρᵣ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )



    end


    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    for face in faces_boundary
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        pₙ = cells[face.owner].var[👉.p]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        ΔS = face.ΔS

        uₙ = 0.0
        vₙ = 0.0
        wₙ = 0.0
        Uₙ = 0.0
        α₁ₙ = cells[face.owner].var[👉.α₁]

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [face.x, face.y, face.z]
        ΔLR = norm(centerᵣ - centerₗ)
        d̂ = 👉.Δt / ρₙ
        
        # continuity
        i += 1
        A_vals[i] += 0.0#( d̂ / ΔLR * ΔS )
        i += 1
        A_vals[i] += 0.0#( 0.5 * face.n̂[1] * ΔS )
        i += 1
        A_vals[i] += 0.0#( 0.5 * face.n̂[2] * ΔS )

        
        # x-momentum
        i += 1
        A_vals[i] += 0.5 * face.n̂[1] * ΔS#( 0.5 * face.n̂[1] * ΔS + ρₙ * d̂ / ΔLR * uₙ * ΔS )
        i += 1
        A_vals[i] += 0.0#( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + 0.5 * ρₙ * Uₙ * ΔS )
        i += 1
        A_vals[i] += 0.0#( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )

        
        # y-momentum
        i += 1
        A_vals[i] += 0.5 * face.n̂[2] * ΔS#( 0.5 * face.n̂[2] * ΔS + ρₙ * d̂ / ΔLR * vₙ * ΔS )
        i += 1
        A_vals[i] += 0.0#( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        i += 1
        A_vals[i] += 0.0#( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + 0.5 * ρₙ * Uₙ * ΔS )



        B[ijStartₗ + 1] -= ( Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )

        
        
    end
 


    A = sparse(A_rows,A_cols,A_vals)
    #x = zeros(Float64, length(cells), 1)

    #println(A)
    #println(B)

    #spy(A, marker=".", markersize=1)
    #gui()

    ps = MKLPardisoSolver()
    ΔQ = solve(ps, A, B)
    


    relax = 1.0

#=
    diagon = 1

    for cell in cells
        istart = (diagon-1)*3
        println(A_vals[istart+1]," ",B[istart+1])
        println(A_vals[istart+2]," ",B[istart+2])
        println(A_vals[istart+3]," ",B[istart+3])

        diagon += 1
    end
=#

    #println("i,j,A")
    #for i in 1:length(A_vals)
        #println(A_rows[i]-1," ",A_cols[i]-1," ",A_vals[i])
    #end
    #println("B")
    #for i in B
        #println(i)
    #end

    #for i in ΔQ

    #    println(i)

    #end

    #sleep(1000.0)




    diagon = 1
    for cell in cells
        
        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        cell.var[👉.p] += relax * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax * ΔQ[ijStart + 3]
        
        diagon += 1
    end


    #return log10(norm(ΔU))
    return log10(norm(ΔQ))
   

end




