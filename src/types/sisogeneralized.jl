ExprLike = Union{Expr,Number,Symbol}

## User should just use TransferFunction
immutable SisoGeneralized <: SisoTf
    expr::ExprLike
    function SisoGeneralized(expr::ExprLike)
        if isa(expr, Expr) && length(expr.args) == 3 && expr.args[1] == :(+) && expr.args[2] == 0
            #Get rid of the zero
            expr = expr.args[3]
        end
        new(expr)
    end
end

SisoGeneralized(str::AbstractString) = SisoGeneralized(parse(str))

function minreal(sys::SisoGeneralized, eps::Real=sqrt(eps()))
    warn("minreal is not implemented for generalized transferfunctions, returning same system")
    sys
end

function print_siso(io::IO, t::SisoGeneralized, var=:s)
    println(io, t.expr)
end

Base.promote_rule{T<:Real}(::Type{SisoGeneralized}, ::Type{T}) = SisoGeneralized
Base.convert(::Type{SisoGeneralized}, b::Real) = SisoGeneralized(b)

Base.zero(::Type{SisoGeneralized}) = SisoGeneralized(0)
Base.zero(::SisoGeneralized) = Base.zero(SisoGeneralized)

Base.length(t::SisoGeneralized) = error("length is not implemented for generalized transferfunctions")
Base.num(t::SisoGeneralized) = error("num is not implemented for generalized transferfunctions")
Base.den(t::SisoGeneralized) = error("den is not implemented for generalized transferfunctions")
pole(t::SisoGeneralized) = error("pole is not implemented for generalized transferfunctions")
tzero(t::SisoGeneralized) = error("tzero is not implemented for generalized transferfunctions")

#This makes sure that the function can compile once
function _preprocess_for_freqresp(sys::SisoGeneralized)
    _f = eval(:(s -> $(sys.expr)))
end

evalfr(f::Function, freq) = f(freq)
evalfr(sys::SisoGeneralized, freq) = _preprocess_for_freqresp(sys)(freq)

function lsimabstract(sys::SisoGeneralized, uin, dt, Tend)
    #TODO make sure U and P are of equal length, fix input arguments, Tend results in half time, make sure u interp is using Tend
    N = round(Int, Tend/dt) + 1
    #T=N*dt
    T = Tend
    dw = pi/T
    omega = linspace(-pi/dt, pi/dt, 2N+1)
    u = [uin; zeros(N)]
    U = fft(u)
    Pf = _preprocess_for_freqresp(sys)
    P = Complex{Float64}[evalfr(Pf, omega[i]*im) for i in 1:2N]
    y = real(ifft(fftshift(P).*U))
    t = dt*(0:N-1)
    y[1:N]
end

function SisoRational(expr::Expr)
    if length(expr.args) == 1
        error("Unexpected operator $(expr.args[1]) in expression")
    end
    if expr.args[1] == :+
        return reduce(+, map(t -> SisoRational(t), expr.args[2:end]))
    elseif expr.args[1] == :*
        return reduce(*, map(t -> SisoRational(t), expr.args[2:end]))
    elseif length(expr.args) == 2 && expr.args[1] == :-
        return - SisoRational(expr.args[2])
    elseif length(expr.args) == 3
        if expr.args[1] == :-
            return SisoRational(expr.args[2]) - SisoRational(expr.args[3])
        elseif expr.args[1] == :/
            return SisoRational(expr.args[2]) / SisoRational(expr.args[3])
        elseif expr.args[1] == :^
            return SisoRational(expr.args[2]) ^ expr.args[3]
        else
            error("Unexpected operator $(expr.args[1]) in expression")
        end
    else
        error("Could not parse \"$(expr.args[1])\" in expression")
    end
end
SisoRational(expr::Symbol) = (expr == :s ? SisoRational([1, 0],[1]) : error("invalid symbol \"$exp\" only \"s\" is allowed"))
SisoRational(expr::Number) = isa(expr,Real) ? SisoRational([expr],[1]) : error("Only real numers are allowed in transferfunction")

==(t1::SisoGeneralized, t2::SisoGeneralized) = (t1.expr == t2.expr)

isapprox(t1::SisoGeneralized, t2::SisoGeneralized; kwargs...) = error("isapprox not implemented for generalized tf")

+(t1::SisoGeneralized, t2::SisoGeneralized) = SisoGeneralized(:($(t1.expr) + $(t2.expr)))
+(t::SisoGeneralized, n::Real) = SisoGeneralized(:($(t.expr) + $n))
+(n::Real, t::SisoGeneralized) = SisoGeneralized(:($n + $(t.expr)))

-(t1::SisoGeneralized, t2::SisoGeneralized) = SisoGeneralized(:($(t1.expr) - $(t2.expr)))
-(n::Real, t::SisoGeneralized) = SisoGeneralized(:($n - $(t.expr)))
-(t::SisoGeneralized, n::Real) = SisoGeneralized(:($(t.expr) - $n))

-(t::SisoGeneralized) = SisoGeneralized(:(- $(t.expr)))

*(t1::SisoGeneralized, t2::SisoGeneralized) = SisoGeneralized(:($(t1.expr) * $(t2.expr)))
*(t::SisoGeneralized, n::Real) = SisoGeneralized(:($(t.expr) * $n))
*(n::Real, t::SisoGeneralized) = SisoGeneralized(:($n * $(t.expr)))

/(n::Real, t::SisoGeneralized) = SisoGeneralized(:($n / $(t.expr)))
/(t::SisoGeneralized, n::Real) = SisoGeneralized(:($(t.expr) / $n))
/(t1::SisoGeneralized, t2::SisoGeneralized) = SisoGeneralized(:($(t1.expr) / $(t2.expr)))
