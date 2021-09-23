
function transport!(
    👉::controls, cells::Vector{mesh.Cell}
)


    for cell in cells

        cell.var[👉.μ] = cell.var[👉.α₁] * 0.001 + cell.var[👉.α₂] * 1.e-5

    end


end


