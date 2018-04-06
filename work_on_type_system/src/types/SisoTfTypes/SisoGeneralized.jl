ExprLike = Union{Expr,Number,Symbol}

# NOTE: Does it make sense to have a type for SisoGeneralized?
## User should just use TransferFunction
immutable SisoGeneralized{T} <: SisoTf{T}
    expr::ExprLike
    function SisoGeneralized{T}(expr::ExprLike) where T
        if isa(expr, Expr) && length(expr.args) == 3 && expr.args[1] == :(+) && expr.args[2] == 0
            #Get rid of the zero
            expr = expr.args[3]
        end
        new{T}(expr)
    end
end

SisoGeneralized(str::AbstractString) = SisoGeneralized(parse(str))
