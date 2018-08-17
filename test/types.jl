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
to_complex(x::Type{DualNumbers.Dual{T}}) where {T} = DualNumbers.Dual{to_complex(typeof(x))}



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
all_fuctions_types = Dict{Function, Tuple{Int, Tuple, String}}(
balance => (   1,  (
    ((BlasFloat,), (T) -> (Matrix{T},), (T) -> (NTuple{3, Matrix{T}})),
                ), "Base.LinAlg.LAPACK.gebal!"),
care => (      1,  (
    ((BlasFloat,), (T) -> NTuple{3, Matrix{T}}, (T) -> Matrix{T}),
                ), "inv, schurfact, ordschur" ),
dare => (      1,  (
    ((BlasFloat,), (T) -> NTuple{3, Matrix{T}}, (T) -> Matrix{T}),
                ), "inv, schurfact, ordschur" ),
dlyap => (     1,  (
    ((BlasFloat,), (T) -> NTuple{2, Matrix{T}}, (T) -> Matrix{T}),
               ), "" ),
lqr => (       1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Matrix{T}, Matrix{T}}, (T) -> Matrix{T}),
                ), "care, dare" ),
dlqr => (      1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
                ), "dare" ),
kalman => (    1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
                ), "care, dare" ),
# TODO remove dkalman?
dkalman => (   1,  (
    ((BlasFloat,), (T) -> NTuple{4, Matrix{T}}, (T) -> Matrix{T}),
                ), "dlqr" ),
# Complicated
# (lqg, ),
# (lqgi, ),
covar => (     1,  (
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Matrix{T}}, (T) -> Matrix{T}),
                ), "Base.lyap"),
norm => (      1,  (
    ((Real,), (T) -> Tuple{CustomLTISystem{T}, Real}, (T) ->  promote_type(Float32, typeof(one(T)/norm(one(T)))) ),
                ), "svdvals, eigvals, type unstable?, unlear if working with complex"),
norminf => (   1,  (
    ((Real,), (T) -> Tuple{CustomLTISystem{T}}, (T) ->  begin
                        T2 = promote_type(Float32, typeof(one(T)/norm(one(T))))
                        return  Tuple{T2,T2}
                    end ),
                ), "svdvals, eigvals, type unstable?, unlear if working with complex"),
gram => (      1,  (
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Symbol}, (T) -> Matrix{T}),
                ), "Base.lyap, dlyap"),
ctrb => (      1,  (
    ((Number,), (T) -> NTuple{2, AbstractMatrix{T}}, (T) -> Matrix{T}),
    ((Number,), (T) -> Tuple{StateSpace{T}}, (T) -> Matrix{T}),
                ), ""),
obsv => (      1,  (
    ((Number,), (T) -> NTuple{2, AbstractMatrix{T}}, (T) -> Matrix{T}),
    ((Number,), (T) -> Tuple{StateSpace{T}}, (T) -> Matrix{T}),
                ), ""),
# TODO Fix acker
place => (     1,  (
    ((RealBlasFloat,), (T) -> Tuple{Matrix{T}, Matrix{T}, Vector{T}}, (T) -> Matrix{promote_type(T,Float64)}),
    ((RealBlasFloat,), (T) -> Tuple{StateSpace{T}, Vector{T}}, (T) -> Matrix{promote_type(T,Float64)}),
                ), "acker, unclear for Complex"),
# Model Simplification
reduce_sys => (1,  (
    ((BlasFloat,), (T) -> Tuple{Matrix{T},Matrix{T},Matrix{T},Matrix{T},BlasFloat}, (T) -> Matrix{T}),
                ), "qrfact"),
sminreal => (  1,  (
    ((Number,), (T) -> Tuple{StateSpace{T}}, (T) -> StateSpace{T}, ""),
                ), ""),
# TODO try to be more specific here, intype=outtype
minreal => (   1,  (
    ((Number,), (T) -> Tuple{CustomLTISystem{T}, Real}, (T) -> CustomLTISystem{T}),
                ), "eigvals (Polynomials.roots)"),
balreal => (   1,  (
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}}, (T) -> Tuple{StateSpace{T}, Matrix{T}}),
                ), "gram, chol, inv"),
baltrunc => (  1,  (
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}}, (T) -> Tuple{StateSpace{T}, Matrix{T}}),
                ), "balreal"),
# # Stability Analysis
# (isstable,  1,  (
#     ((,),)
#                 ), "pole"),
# TODO figure out appropriate restriction for eigvals (in StateSpace case)
# TODO return type is TR in SisoZpk{T,TR} case, complex(T) not working for some types
# TODO define complex(DualNumbers.Dual{T<:BlasFloat}) = DualNumbers.Dual{Complex{T}}
pole => (      1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Vector{complex(T)}),
    ((BlasFloat,), (T) -> TransferFunction{SisoRational{T}}, (T) -> Vector{complex(T)}),
    ((Number,), (T) -> TransferFunction{SisoZpk{T,<:Number}}, (T) -> Vector{complex(T)}),
                ), "eigvals (Polynomials.roots)"),
# TODO Not correct for Zpk, may have real zeros/poles
tzero => (     1,  (
    ((BlasFloat,), (T) -> StateSpace{T}, (T) -> Vector{complex(T)}),
    ((BlasFloat,), (T) -> TransferFunction{SisoRational{T}}, (T) -> Vector{complex(T)}),
    ((Number,), (T) -> TransferFunction{SisoZpk{T,<:Number}}, (T) -> Vector{to_complex(T)}),
                ), "Base.LinAlg.LAPACK.gebal!"),
# TODO include other input arguments, some one(T) should be typeof(s)
dcgain => (    1,  (
    ((Number,), (T) -> CustomLTISystem{T}, (T) -> Base.promote_op(/,T,T)),
                ), "evalfr"),
#complicated
#(zpkdata, ),
# TODO not really true for Zpk, and too general for SS, SisoRational
damp => (      1,  (
    ((Number,), (T) -> LTISystem, (T) -> NTuple{3,Vector{float(T)}}),
                ), "cos, angle, log, poles"),
#(dampreport, ),
markovparam => (1, (
    ((Number,), (T) -> StateSpace{T}, (T) -> Matrix{T},)
                ), "" ),

# TODO Margin issue 146
#(margin, ),
#(delaymargin, "margin"),
#(gangoffour, ),
# Connections
append => (    2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T1},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type(T1,T2)}} ),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{promote_type(T1,T2)}}),
                ), ""),
series => (    2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T1},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type(T1,T2)}} ),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{promote_type(T1,T2)}}),
                ), ""),
parallel => (  2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T1},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type(T1,T2)}} ),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{promote_type(T1,T2)}}),
                ), ""),
# No support for DualNumbers with SisoZpk, only real/complex
feedback => ( 2,  (
        ((Number, Number,), (T1,T2) -> Tuple{StateSpace{T1},StateSpace{T2}}, (T1,T2) -> StateSpace{promote_type_div(T1,T2)}),
        ((Number, Number,), (T1,T2) -> Tuple{TransferFunction{SisoRational{T1}},TransferFunction{SisoRational{T2}}},
            (T1,T2) ->  TransferFunction{SisoRational{promote_type_div(T1,T2)}} ),
        ((Union{Real, Complex}, Union{Real, Complex},), (T1,T2) -> Tuple{TransferFunction{SisoZpk{T1}},TransferFunction{SisoZpk{T2}}},
            (T1,T2) ->  TransferFunction{SisoZpk{complex(promote_type_div(T1,T2))}}),
                ), "Polynomials.roots"),
feedback2dof => (1,(
        ((Number,), (T) -> Tuple{TransferFunction{SisoRational{T}}, Vector{T}, Vector{T}, Vector{T}},
            (T) ->  TransferFunction{SisoRational{T}} ),
        ((Union{Real, Complex},), (T) -> Tuple{TransferFunction{SisoZpk{T}}, Vector{T}, Vector{T}, Vector{T}},
            (T) ->  TransferFunction{SisoZpk{T, complex(T)}} ),
        ), "Polynomials.roots" ),
# Discrete
c2d => (   1,  (
        ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Number, Symbol}, (T) -> StateSpace{T}),
        ), "expm"),
# Time Response
# discrete/continuous make type-unstable because of c2d
# Only correct for discrete below ?
# TODO test with function and :method options
# Output: y,t,x
step => (  1,  (
        ((Number,), (T) -> Tuple{StateSpace{T}, Vector{T}},
            T -> Tuple{Array{Float64,3}, Vector{T}, Array{Float64,3}}),
        ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{T}},
            T -> Tuple{Array{Float64,3}, Vector{T}, Array{Float64,3}}),
        ), "c2d, "),
impulse => (1,  (
        ((Number,), (T) -> Tuple{StateSpace{T}, Vector{T}},
            T -> Tuple{Array{Float64,3}, Vector{T}, Array{Float64,3}}),
        ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{T}},
            T -> Tuple{Array{Float64,3}, Vector{T}, Array{Float64,3}}),
        ), "c2d, "),
lsim => (  1,  (
        ((Number,), (T) -> Tuple{StateSpace{T}, Matrix{T}, Vector{T}},
            T -> Tuple{Matrix{Float64},Array{T},Matrix{Float64}}),
        ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Matrix{T}, Vector{T}},
            T -> Tuple{Matrix{Float64}, Matrix{T}, Matrix{Float64}}),
        ), "c2d, "),

# TODO solve and Simulator have no type information
# (solve, ),
# (Simulator, ),

# Frequency Response
# Not completely accurate, if w has higher accuracy, that will be reflected in output
freqresp => (1,  (
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}}),
        (T) -> Array{promote_type(T, Complex64),3}),
    ((Number,), (T) -> (Tuple{StateSpace{T}, T}),
        (T) -> Array{promote_type(T, Complex64),3}),
        ), "evalfr, (exp for discrete)"),
# Roughly true

evalfr => (1,  (
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Real}),
        (T) -> Base.promote_op(/, T, T)),
    ((Number,), (T) -> (Tuple{StateSpace{T}, Real}),
        (T) -> Base.promote_op(/, T, T)),
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Complex}),
        (T) -> promote_type_div(T, complex(T))),
    ((Number,), (T) -> (Tuple{StateSpace{T}, Complex}),
        (T) -> promote_type_div(T, complex(T))),
        ), ""),

bode => (  1,  (
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}}),
        (T) -> Tuple{Array{promote_type(real(T), Float32),3},
                     Array{promote_type(real(T), Float32),3},
                        Vector{<:Real}}),
    ((Number,), (T) -> (Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}}),
        (T) -> Tuple{Array{promote_type(real(T), Float32),3},
                     Array{promote_type(real(T), Float32),3},
                        Vector{<:Real}}),
        ), "freqresp"),
nyquist => (  1,  (
    ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}},
        (T) -> Tuple{Array{promote_type(real(T), Float32),3},
                     Array{promote_type(real(T), Float32),3},
                        Vector{<:Real}}),
    ((Number,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}},
        (T) -> Tuple{Array{promote_type(real(T), Float32),3},
                     Array{promote_type(real(T), Float32),3},
                        Vector{<:Real}}),
        ), "freqresp"),
# TODO only Float64 output!
sigma => ( 1,  (
    ((BlasFloat,), (T) -> Tuple{TransferFunction{<:SisoTf{T}}, Vector{<:Real}},
        (T) -> Tuple{Array{Float64,2},Vector{<:Real}}),
    ((BlasFloat,), (T) -> Tuple{StateSpace{T}, Vector{<:Real}},
        (T) -> Tuple{Array{Float64,2},Vector{<:Real}}),
        ), "freqresp, svdvals"),

# # utilities
# (numpoly, ),
# (denpoly, ),
)


ss1matrix(T) = (T[-1 1; 0 1], T[1 0;0 1], T[1 0], fill(T(0),1,2))
ss2matrix(T) = (T[-3 1; 0 1], T[1 0;0 1], T[1 0], fill(T(0),1,2))
ss3matrix(T) = (T[-1 1; 0 -1], T[1 0;0 1], T[1 0], fill(T(0),1,2))
ss4matrix(T) = (T[-1 1; 0 -1], reshape(T[0,1],2,1), T[1 0], fill(T(0),1,1))

tf1matrix(T) = [SisoRational{T}([1,-1],[1,0,-1]) SisoRational{T}([1,],[1,0,-1])]
tf2matrix(T) = [SisoRational{T}([1,-1],[1,2,-3]) SisoRational{T}([1,],[1,2,-3])]

zpk1matrix(T) = [SisoZpk{T,Complex{T}}([-1,],[1, -1],1) SisoZpk{T,Complex{T}}([],[1,-1],1)]
zpk2matrix(T) = [SisoZpk{T,Complex{T}}([-1,],[-3, 1],1) SisoZpk{T,Complex{T}}([],[-3,1],1)]

systemsdict(T,dt=0) = Dict(
    "statespace1" => StateSpace{T, Matrix{T}}(ss1matrix(T)..., dt),
    "statespace2" => StateSpace{T, Matrix{T}}(ss2matrix(T)..., dt),
    "statespace3" => StateSpace{T, Matrix{T}}(ss3matrix(T)..., dt), # stable mimo
    "statespace4" => StateSpace{T, Matrix{T}}(ss4matrix(T)..., dt), # stable siso
    "tf1" => TransferFunction{SisoRational{T}}(tf1matrix(T), dt),
    "tf2" => TransferFunction{SisoRational{T}}(tf2matrix(T), dt),
    "zpk1" => TransferFunction{SisoZpk{T,Complex{T}}}(zpk1matrix(T), dt),
    "zpk2" => TransferFunction{SisoZpk{T,Complex{T}}}(zpk2matrix(T), dt),
    )

# (ss1 +0im, ss2 +im*ss1,
#  tf1 +0im, tf2 +im*tf2,
#  zpk1+0im, zpk1+im*zpk2 )
systemsdict_complex(T,dt=0) = Dict(
        "statespace1" => StateSpace{Complex{T}, Matrix{Complex{T}}}(
            (ss1matrix(T) .+ 0 .* im.*ss1matrix(T))..., dt),
        "statespace2" => StateSpace{Complex{T}, Matrix{Complex{T}}}(
            (ss2matrix(T) .+ 1 .* im.*ss1matrix(T))..., dt),
        "tf1" => TransferFunction{SisoRational{Complex{T}}}(
            (tf1matrix(T) .+ 0 .* im.*tf1matrix(T)), dt),
        "tf2" => TransferFunction{SisoRational{Complex{T}}}(
            (tf2matrix(T) .+ 1 .* im.*tf1matrix(T)), dt),
        #"zpk1" => TransferFunction{SisoZpk{T,Complex{T}}}(
        #    (zpk1matrix(T) .+ 0 .* im.*zpk1matrix(T))..., 0),
        #"zpk2" => TransferFunction{SisoZpk{T,Complex{T}}}(
        #    (zpk2matrix(T) .+ 1 .* im.*zpk1matrix(T))..., 0),
    )

sysoftypes = [
    (Float32, systemsdict(Float32,0)),
    (Float64, systemsdict(Float64,0)),
    (Int, systemsdict(Int,0)),
    (Complex{Float32}, systemsdict_complex(Float32,0)),
    (Complex{Float64}, systemsdict_complex(Float64,0)),
]

valsdict = Dict(
    "tvec1" => collect(1:100),
    "freqvec1" => logspace(-1,1,50),
    "covarmat1" => [1 2; 2 1] ,
)



# Not able to infer type-stability without this
function infer_type(fun, args, kwargs...)
    dummy_call(fun, args, kwargs...) = fun(args...; kwargs...)
    return Base.return_types(dummy_call, typeof((fun, args, kwargs...)))[1]
    #return Base.code_typed(dummy_call, typeof((fun, args, kwargs...)))[1][2]
end

timevec1(T) = (real(T).(valsdict["tvec1"]),)
covarmat1(T) = (real(T).(valsdict["covarmat1"]))

functions_and_args = [
    (step,      timevec1,       Dict(:method=>:zoh), nothing),
    (freqresp,  timevec1,       Dict(), nothing),
    (bode,      timevec1,       Dict(), nothing),
    (nyquist,   timevec1,       Dict(), nothing),
    (sigma,     timevec1,       Dict(), nothing),
    (covar,     covarmat1,      Dict(), nothing),
    (norm,      (T) -> (2,) ,   Dict(), nothing),
    (norm,      (T) -> (Inf,),  Dict(), nothing),
    (norminf,   (T) -> (),      Dict(), nothing),
    (gram,      (T) -> (:c,),   Dict(), isstable),
    (gram,      (T) -> (:o,),   Dict(), isstable),
    (ctrb,      (T) -> (),      Dict(), nothing),
    (obsv,      (T) -> (),      Dict(), nothing)
]

function run_test(sysoftypes, all_fuctions_types, functions_and_args)
    @testset "Testing type stability" begin
        @testset "for function: $fun" for
            (fun, argfun, kwargs, requrement) in functions_and_args
            @testset "with T=$T" for (T, sysdict) in sysoftypes
                fun_defs = all_fuctions_types[fun]
                found_fun_def = false
                for (Tset,inT,outT) in fun_defs[2]
                    if Tuple{(T,)...} <: Tuple{Tset...} # Check that this T is allowed for function
                        found_fun_def = true
                        for (key,sys) in sysdict        # Get a candidate system
                            argin = (sys, argfun(T)...)
                            if requrement == nothing || requrement(sys) # Test if function supports system
                                if isa(argin, inT(T))       # Check if we support this system with argin
                                    @testset "with sys=$key of tyope?$(typeof(sys))" begin
                                        out = fun(argin...; kwargs...)
                                        out_T_actual = typeof(out)
                                        @test out_T_actual <: outT(T)
                                        if ! (out_T_actual <: outT(T))  # Check output against expected output type
                                            println("In fun=$fun, T=$T, Tset=$Tset, syskey=$key, typeof(sys) = $(typeof(sys))\ngot\n$(out_T_actual)\nbut expected\n$(outT(T))")
                                        end
                                        inferT = infer_type(fun, argin, kwargs...)
                                        @test isleaftype(inferT)
                                        if !isleaftype(inferT)          # Check type stability
                                            println("In fun=$fun, T=$T, Tset=$Tset, typeof(sys) = $(typeof(sys))\n Type instability, got\n $inferT")
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                if !found_fun_def
                    println("No definition found for $fun with T=$T")
                end
            end
        end
    end
end
