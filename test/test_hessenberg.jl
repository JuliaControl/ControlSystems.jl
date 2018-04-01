@testset "test_hessenberg" begin

srand(31415)

A = randn(100,100)
F = hessfact(A)
w = 2.0
B = randn(100)

res1 = (I*w*im - A)\B
res2 = solvehess(F, B, w)

@test res1 ≈ res2 rtol = 1e-12

A = complex(randn(100,100))
F = hessfact(A)
w = 2
B = randn(100)

res1 = (I*w*im - A)\B
res2 = solvehess(F, B, w)

@test res1 ≈ res2 rtol = 1e-12

A = randn(10,10)
B = randn(10,2)
C = randn(3,10)
D = randn(3,2)

w = 0.2
sys1 = ss(A,B,C,D)
sys2 = ControlSystems._preprocess_for_freqresp(sys1)::HessenbergSS

res0 = D + C*((I*im*w-A)\B)
res1 = evalfr(sys1,im*w)
res2 = evalfr(sys2,im*w)

@test res0 ≈ res1 rtol = 1e-12
@test res1 ≈ res2 rtol = 1e-12

end
