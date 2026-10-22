using ControlSystemsBase, LinearAlgebra, Printf

# Modal / pole-residue evaluation (Float64), returns response + cond(V)
function modal_eval(A,B,C,D,w)
    E = eigen(A)
    λ = ComplexF64.(E.values); V = ComplexF64.(E.vectors)
    ctil = vec(complex.(C)*V); btil = V \ ComplexF64.(B[:,1])
    res = ctil .* btil
    R = ComplexF64[ D[1,1] + sum(res ./ (complex(0.0,ω) .- λ)) for ω in w ]
    R, cond(V)
end

# High-precision reference: resolvent in Complex{BigFloat}
function ref_eval(A,B,C,D,w)
    setprecision(256) do
        Ab=Complex{BigFloat}.(A); Bb=Complex{BigFloat}.(B[:,1])
        Cb=Complex{BigFloat}.(C); d=Complex{BigFloat}(D[1,1])
        [ (Cb*((complex(BigFloat(0),BigFloat(ω))*I - Ab)\Bb))[1] + d for ω in w ]
    end
end

relerr(R, Rref) = maximum(abs.(ComplexF64.(R) .- ComplexF64.(Rref))) / maximum(abs, ComplexF64.(Rref))

function report(name, sys, w)
    A,B,C,D = sys.A, sys.B, sys.C, sys.D
    Rref = ref_eval(A,B,C,D,w)
    Rh   = vec(freqresp(sys, w))                 # package Hessenberg path
    Rm, κ = modal_eval(A,B,C,D,w)
    @printf("%-32s nx=%2d  cond(V)=%.2e   Hessenberg relerr=%.2e   modal relerr=%.2e\n",
            name, sys.nx, κ, relerr(Rh,Rref), relerr(Rm,Rref))
end

w = exp10.(LinRange(-2, 2, 60))

report("random (normal-ish)", ssrand(1,1,20), w)
report("8 repeated poles at -1", ss(tf(1.0,[1.0,1.0]))^8, w)
report("12 repeated poles at -1", ss(tf(1.0,[1.0,1.0]))^12, w)
# clustered (non-defective but near-defective) eigenvalues via upper-triangular A
let n=20
    A = triu(0.1 .* randn(n,n)); A[diagind(A)] .= -1.0 .- 0.001 .* (1:n)  # tightly clustered evals
    B = randn(n,1); C = randn(1,n); D = zeros(1,1)
    report("clustered eigenvalues (triu)", ss(A,B,C,D), w)
end
# lightly damped pole near imaginary axis, evaluate across it
report("lightly damped (zeta=1e-4)", ss(tf(1.0,[1.0, 2e-4, 1.0])), exp10.(LinRange(-1,1,200)))
