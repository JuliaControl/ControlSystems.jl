RealBlasFloat = Union{Float32, Float64}
ComplexBlasFloat = Union{Complex{Float32}, Complex{Float64}}
import Base.LinAlg: BlasFloat
# Constructors
LTISystem,
StateSpace,
TransferFunction,
ss,
tf,
zpk,
ss2tf,
LQG,
primitivetype,

# Linear Algebra
# (function, (nin, nout), NTuple{Tuple{TS, TFin, TFout}, note)
# where TS is type of allowable inputs, e.g Tuple{Matrix{T}, T} where T<:BalsFloat
# TFin, function from TS to input types
# TFout, function ftom TS to output types
# Note, String explaning bottleneck or other problems
(
(balance,   1,  (
    ((BlasFloat,), (T) -> (Matrix{T},), (T) -> (NTuple{3, Matrix{T}})),
                ), "Base.LinAlg.LAPACK.gebal!"),
(care,      1,  (
    ((BlasFloat,), (T) -> NTuple{3, Matrix{T}}, (T) -> Matrix{T}),
                ), "inv, schurfact, ordschur" ),
(dare,      1,  (
    ((BlasFloat,), (T) -> NTuple{3, Matrix{T}}, (T) -> Matrix{T}),
                ), "inv, schurfact, ordschur" ),
(dlyap,     1,  (
    ((BlasFloat,), (T) -> NTuple{2, Matrix{T}}, (T) -> Matrix{T}),
               ), "" ),
(lqr,       1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Matrix{T}, Matrix{T}}, (T) -> Matrix{T}),
                ), "care, dare" ),
(dlqr,      1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
                ), "dare" ),
(kalman,    1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
                ), "care, dare" ),
# TODO remove dkalman?
(dkalman,   1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
                ), "dlqr" ),
# Complicated
# (lqg, ),
# (lqgi, ),
(covar,     1,  (
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Matrix{T}}, (T) -> Matrix{T}),
                ), "Base.lyap"),
(norm,      1,  (
    ((Real,), (T) -> Tuple{LTISystem{T}, Real}, (T) ->  begin
                        T2 = promote_type(Float32, typeof(one(T)/norm(one(T))))
                        return  Tuple{T,T}
                    end ),
                ), "svdvals, eigvals, type unstable?, unlear if working with complex"),
(norminf,   1,  (
    ((Real,), (T) -> Tuple{LTISystem{T}, Real}, (T) ->  begin
                        T2 = promote_type(Float32, typeof(one(T)/norm(one(T))))
                        return  Tuple{T,T}
                    end ),
                ), "svdvals, eigvals, type unstable?, unlear if working with complex"),
(gram,      1,  (
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Symbol}, (T) -> Matrix{T}),
                ), "Base.lyap, dlyap"),
(ctrb,      1,  (
    ((Number,), (T) -> NTuple{2, AbstractMatrix{T}}, (T) -> Matrix{T}),
    ((Number,), (T) -> StateSpace{T}, (T) -> Matrix{T}),
                ), ""),
(obsv,      1,  (
    ((Number,), (T) -> NTuple{2, AbstractMatrix{T}}, (T) -> Matrix{T}),
    ((Number,), (T) -> StateSpace{T}, (T) -> Matrix{T}),
                ), ""),
# TODO Fix acker
(place,     1,  (
    ((RealBlasFloat,), (T) -> Tuple{Matrix{T}, Matrix{T}, Vector{T}}, (T) -> Matrix{promote_type(T,Float64)}),
    ((RealBlasFloat,), (T) -> Tuple{StateSpace{T}, Vector{T}}, (T) -> Matrix{promote_type(T,Float64)}),
                ), "acker, unclear for Complex"),
# Model Simplification
(reduce_sys,1,  (
    ((BlasFloat,), (T) -> Tuple{Matrix{T},Matrix{T},Matrix{T},Matrix{T},BlasFloat}, (T) -> Matrix{T}),
                ), "qrfact"),
(sminreal,  1,  (
    ((Number,), (T) -> StateSpace{T}, (T) -> StateSpace{T}, ""),
                ), ""),
# TODO try to be more specific here, intype=outtype
(minreal,   1,  (
    ((Number,), (T) -> Tuple{LTISystem{T}, Real}, (T) -> LTISystem{T}),
                ), "eigvals (Polynomials.roots)"),
(balreal,   1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Tuple{StateSpace{T}, Matrix{T}}),
                ), "gram, chol, inv"),
(baltrunc,  1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Tuple{StateSpace{T}, Matrix{T}}),
                ), "balreal"),
# Stability Analysis
(isstable,  1,  (
    ((,),)
                ), "pole"),
# TODO figure out appropriate restriction for eigvals (in StateSpace case)
# TODO return type is TR in SisoZpk{T,TR} case, complex(T) not working for some types
# TODO define complex(DualNumbers.Dual{T<:BlasFloat}) = DualNumbers.Dual{Complex{T}}
(pole,      1, (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Vector{complex(T)}),
    ((BlasFloat,), (T) -> TransferFunction{SisoRational{T}}, (T) -> Vector{complex(T)}),
    ((Number,), (T) -> TransferFunction{SisoZpk{T,<:Number}}, (T) -> Vector{complex(T)}),
                ), "eigvals (Polynomials.roots)"),
(tzero,     1, (
    (())
)),
(dcgain, ),
(zpkdata, ),
(damp, ),
(dampreport, ),
(markovparam, ),
(margin, ),
(delaymargin, ),
(gangoffour, ),
# Connections
(append, ),
(series, ),
(parallel, ),
(feedback, ),
(feedback2dof, ),
# Discrete
(c2d, ),

# Time Response
(step, ),
(impulse, ),
(lsim, ),
(solve, ),
(Simulator, ),

# Frequency Response
(freqresp, ),
(evalfr, ),
(bode, ),
(nyquist, ),
(sigma, ),

# utilities
(numpoly, ),
(denpoly, ),
