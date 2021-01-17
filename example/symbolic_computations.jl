using Revise
using ControlSystems
using SymPy

# Some basic demonstations of working with symbolic LTI systems
# Functionality is rather limited, and for complicated expressions the
# printing is awful.
# For additional functionality, see https://github.com/JuliaControl/SymbolicControlSystems.jl


# Need to modify some functions, to work with the Sym type
# Some of these changes are on their way into the respective packages
Base.complex(::Type{SymPy.Sym}) = Sym
Base.eigvals(a::Matrix{Sym}) = Vector{Sym}(collect(keys(call_matrix_meth(a, :eigenvals))))

function Polynomials.roots(p::Polynomials.Polynomial{SymPy.Sym})
    x = Sym("x")
    return SymPy.solve(p(x), x)
end

function SymPy.simplify(G::TransferFunction{ControlSystems.SisoZpk{Sym,Sym}})
    z, p, k = zpkdata(G)

    return zpk(map(x->simplify.(x), z), map(x->simplify.(x), p), map(x->simplify.(x), k))
end

function ControlSystems.impulse(sys::StateSpace{Sym})
    t = symbols("t", real=true)
    return simplify.(sys.C * exp(sys.A*t) * sys.B)
end


# Define a symbolic parameter
a = symbols("a", real=true)

# Define a statespace and a trasnfer function
sys = ss([1 a; a 1], [0; 1], [1 0], 0)

s = tf("s")
G = (s/(5a)) / ((s/a)^2 + s/(5a) + 1)


# Simple conversions
@edit zpk(sys)
tf(sys)

# The coefficient looks a bit complicated, but simplifying gives..
z, p, k = zpkdata(tf(sys))
simplify.(k[1])

zpkdata(G)
ss(G)


# Controllability/observability matrices
Mc = ctrb(sys)
Mo = obsv(sys)


## Compute the L∞ norm (or actually  L∞ norm squared) for two symbolic systems
w = symbols("w", real=true)

sys_to_consider = sys # G

sys_fr = simplify(evalfr(sys_to_consider, im*w)[1])
sys_fr_mag = simplify(abs(sys_fr)^2)

n, _ = fraction( diff(sys_fr_mag, w) )
roots = SymPy.solve(SymPy.Poly(n), w)

real_roots = roots[SymPy.imag(roots) .== 0]

maximum([subs(sys_fr_mag, w => r) for r in real_roots])


# Compute the impulse resonse of some systems (on statespace form)
impulse(sys)[1]
simplify(impulse(ss(G))[1])

rosenbrock = [1/(s+1) 1/(s+a); 1/(s+1) 1/(s+1)]
ss(rosenbrock)
impulse(ss(rosenbrock))


# Analytic impulse response
sys2 = ss([-2.5 0;1 1.5],[1;3],[1 2],Sym(2.5))
impulse(sys2)[1]
