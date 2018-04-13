function tfg(systems::Array, Ts::Real=0; kwargs...)
    ny, nu = size(systems, 1, 2)

    T = Float64
    matrix = Matrix{SisoGeneralized{T}}(ny, nu)
    for o=1:ny
        for i=1:nu
            matrix[o, i] = SisoGeneralized{T}(systems[o, i])
        end
    end
    return TransferFunction(matrix, Float64(Ts))
end

tfg(var::Union{AbstractString,ExprLike}, Ts=0; kwargs...) = tfg([var], Ts; kwargs...)
