using ControlSystems, DualNumbers

RealBlasFloat = Union{Float32, Float64}
ComplexBlasFloat = Union{Complex{Float32}, Complex{Float64}}
import Base.LinAlg: BlasFloat
import ControlSystems: SisoTf, SisoRational, SisoZpk
# TODO maybe define LTISystem like this?
const CustomLTISystem{T} = Union{TransferFunction{<:SisoTf{T}}, StateSpace{T}}
function promote_type_div(T1,T2)
    T0 = promote_type(T1,T2)
    return Base.promote_op(/, T0, T0)
end

to_complex(x::Type{<:Number}) = complex(typeof(x))
to_complex(x::Type{DualNumbers.Dual{T}}) where T = DualNumbers.Dual{to_complex(typeof(x))}



# # Constructors
# LTISystem,
# StateSpace,
# TransferFunction,
# ss,
# tf,
# zpk,
# ss2tf,
# LQG,
# primitivetype,

# Linear Algebra
# (function, (nin, nout), NTuple{Tuple{TS, TFin, TFout}, note)
# where TS is type of allowable inputs, e.g Tuple{Matrix{T}, T} where T<:BalsFloat
# TFin, function from TS to input types
# TFout, function ftom TS to output types
# Note, String explaning bottleneck or other problems
all_fuctions_types = [
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
    ((Real,), (T) -> Tuple{CustomLTISystem{T}, Real}, (T) ->  begin
                        T2 = promote_type(Float32, typeof(one(T)/norm(one(T))))
                        return  Tuple{T,T}
                    end ),
                ), "svdvals, eigvals, type unstable?, unlear if working with complex"),
(norminf,   1,  (
    ((Real,), (T) -> Tuple{CustomLTISystem{T}, Real}, (T) ->  begin
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
    ((Number,), (T) -> Tuple{CustomLTISystem{T}, Real}, (T) -> CustomLTISystem{T}),
                ), "eigvals (Polynomials.roots)"),
(balreal,   1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Tuple{StateSpace{T}, Matrix{T}}),
                ), "gram, chol, inv"),
(baltrunc,  1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Tuple{StateSpace{T}, Matrix{T}}),
                ), "balreal"),
# # Stability Analysis
# (isstable,  1,  (
#     ((,),)
#                 ), "pole"),
# TODO figure out appropriate restriction for eigvals (in StateSpace case)
# TODO return type is TR in SisoZpk{T,TR} case, complex(T) not working for some types
# TODO define complex(DualNumbers.Dual{T<:BlasFloat}) = DualNumbers.Dual{Complex{T}}
(pole,      1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Vector{complex(T)}),
    ((BlasFloat,), (T) -> TransferFunction{SisoRational{T}}, (T) -> Vector{complex(T)}),
    ((Number,), (T) -> TransferFunction{SisoZpk{T,<:Number}}, (T) -> Vector{complex(T)}),
                ), "eigvals (Polynomials.roots)"),
# TODO Not correct for Zpk, may have real zeros/poles
(tzero,     1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Vector{complex(T)}),
    ((BlasFloat,), (T) -> TransferFunction{SisoRational{T}}, (T) -> Vector{complex(T)}),
    ((Number,), (T) -> TransferFunction{SisoZpk{T,<:Number}}, (T) -> Vector{to_complex(T)}),
                ), "Base.LinAlg.LAPACK.gebal!"),
# TODO include other input arguments, some one(T) should be typeof(s)
(dcgain,    1,  (
    ((Number,), (T) -> CustomLTISystem{T}, (T) -> Base.promote_op(/,T,T)),
                ), "evalfr"),
#complicated
#(zpkdata, ),
# TODO not really true for Zpk, and too general for SS, SisoRational
(damp,      1,  (
    ((Number,), (T) -> LTISystem, (T) -> NTuple{3,Vector{float(T)}}),
                ), "cos, angle, log, poles"),
#(dampreport, ),
(markovparam,1, (
    ((Number,), (T) -> StateSpace{T}, (T) -> Matrix{T},)
                ), "" ),

# TODO Margin issue 146
#(margin, ),
#(delaymargin, "margin"),
#(gangoffour, ),
# Connections
(append,    2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T2},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type(T1,T2)}} ),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{promote_type(T1,T2)}}),
                ), ""),
(series,    2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T2},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type(T1,T2)}} ),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{promote_type(T1,T2)}}),
                ), ""),
(parallel,  2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T2},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type(T1,T2)}} ),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{promote_type(T1,T2)}}),
                ), ""),
# No support for DualNumbers with SisoZpk, only real/complex
(feedback, 2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T2},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type_div(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type_div(T1,T2)}} ),
        ((Union{Real, Complex}, Union{Real, Complex},), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{complex(promote_type_div(T1,T2))}}),
                ), "Polynomials.roots"),
(feedback2dof,1,(
        ((Number,), (T) -> Tuple{TransferFunction{SisoRational{T}}, Vector{T}, Vector{T}, Vector{T}},
            (T) ->  TransferFunction{SisoRational{T}} ),
        ((Union{Real, Complex},), (T) -> Tuple{TransferFunction{SisoZpk{T}}, Vector{T}, Vector{T}, Vector{T}},
            (T) ->  TransferFunction{SisoZpk{T, complex(T)}} ),
        ), "Polynomials.roots" ),
# Discrete
(c2d,   1,  (
        ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Number, Symbol}, (T) -> StateSpace{T}),
        ), "expm"),
# Time Response
# discrete/continuous make type-unstable because of c2d
# Only correct for discrete below ?
# TODO test with function and :method options
# Output: y,t,x
(step,  1,  (
        ((Number,), (T) -> Tuple{StateSpace{T}, Vector{T}},
            T -> Tuple{Matrix{T}, Vector{T}, Matrix{T}}),
        ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{T}},
            T -> Tuple{Matrix{Base.promote_op(/,T,T)}, Vector{T}, Matrix{Base.promote_op(/,T,T)}}),
        ), "c2d, "),
(impulse,1,  (
        ((Number,), (T) -> Tuple{StateSpace{T}, Vector{T}},
            T -> Tuple{Matrix{T}, Vector{T}, Matrix{T}}),
        ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{T}},
            T -> Tuple{Matrix{Base.promote_op(/,T,T)}, Vector{T}, Matrix{Base.promote_op(/,T,T)}}),
        ), "c2d, "),
(lsim,  1,  (
        ((Number,), (T) -> Tuple{StateSpace{T}, Array{T}, Vector{T}},
            T -> Array{T, 3}),
        ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Array{T}, Vector{T}},
            T -> Array{Base.promote_op(/, T, T), 3}),
        ), "c2d, "),

# TODO solve and Simulator have no type information
# (solve, ),
# (Simulator, ),

# Frequency Response
# Not completely accurate, if w has higher accuracy, that will be reflected in output
(freqresp,1,  (
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}}),
        (T) -> Array{promote_type(T, Complex64),3}),
    ((Number,), (T) -> (Tuple{StateSpace{T}, T}),
        (T) -> Array{promote_type(T, Complex64),3}),
        ), "evalfr, (exp for discrete)"),
# Roughly true

(evalfr,1,  (
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Real}),
        (T) -> Base.promote_op(/, T, T)),
    ((Number,), (T) -> (Tuple{StateSpace{T}, Real}),
        (T) -> Base.promote_op(/, T, T)),
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Complex}),
        (T) -> promote_type_div(T, complex(T))),
    ((Number,), (T) -> (Tuple{StateSpace{T}, Complex}),
        (T) -> promote_type_div(T, complex(T))),
        ), ""),

(bode,  1,  (
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}}),
        (T) -> Tuple{Array{promote_type(T, Float64),3},
                     Array{promote_type(T, Float64),3},
                        Vector{<:Real}}),
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}}),
        (T) -> Tuple{Array{promote_type(T, Float64),3},
                     Array{promote_type(T, Float64),3},
                        Vector{<:Real}}),
        ), "freqresp"),
(nyquist,  1,  (
    ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}},
        (T) -> Tuple{Array{promote_type(T, Float64),3},
                     Array{promote_type(T, Float64),3},
                        Vector{<:Real}}),
    ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}},
        (T) -> Tuple{Array{promote_type(T, Float64),3},
                     Array{promote_type(T, Float64),3},
                        Vector{<:Real}}),
        ), "freqresp"),
# TODO only Float64 output!
(sigma, 1,  (
    ((BlasFloat,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}},
        (T) -> Tuple{Array{Float64,2},Vector{<:Real}}),
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Vector{<:Real}},
        (T) -> Tuple{Array{Float64,2},Vector{<:Real}}),
        ), "freqresp, svdvals"),

# # utilities
# (numpoly, ),
# (denpoly, ),
]


ss1matrix(T) = (T[-1 1; 0 1], T[1 0;0 1], T[1 0], fill(T(0),1,2))
args = [
    (Float32, Dict(
        "statespace" => StateSpace{Float32, Matrix{Float32}}(ss1matrix(Float32)..., 0),
        "tf" => ,
        "zpk" =>,
        )

]
# NOTE
# Test isleaftype(Base.code_typed(f, typeof(args))[1][2])
# for type stability of input -> output
