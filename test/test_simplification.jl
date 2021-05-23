@testset "test_simplification" begin
## SMINREAL ##
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])
@test sminreal(G) == G
@test sminreal(G[1, 1]) == ss(0)
@test sminreal(G[1, 2]) == ss([-1 -2.5; 4 0], [1; 0], [3 0], [0])
@test sminreal(G[2, 1]) == ss([-5], [2], [-2], [1])
@test sminreal(G[2, 2]) == ss([-6], [2], [1], [0])


## MINREAL ##

s = tf("s")
P = [1/(s+1) 2/(s+3); 1/(s+1) 1/(s+1)]
sys = ss(P)
sysmin = minreal(sys)

@test size(sysmin.A,1) == 3 # Test that the reduction of sys worked

@test hinfnorm(sys - sysmin)[1] < 1e-15 # And that the answer is correct

@test_broken balreal(sys-sysmin)

@test all(sigma(sys-sysmin, [0.0, 1.0, 2.0])[1] .< 1e-15)  # Previously crashed because of zero dimensions in tzero

t = 0:0.1:10
y1,x1 = step(sys,t)[[1,3]]
y2,x2 = step(sysmin,t)[[1,3]]
@test sum(abs2,y1.-y2) < 1e-6 # Test that the output from the two systems are the same



import LinearAlgebra.eigsortby
polecompare(sys1, sys2) = isapprox(sort(pole(sys1), by=eigsortby), sort(pole(sys2), by=eigsortby), atol=1e-3)

almost_diagonal(X) = norm(X-diagm(diag(X))) < sqrt(eps())

##
Random.seed!(0)
sys = ssrand(1,1,2)
sysr = minreal(sys)
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

sys = [sys sys] # unobservable
sysr = minreal(sys)
@test sysr.nx == 2
@test polecompare(sys, [sysr sysr])
@test almost_diagonal(gram(sysr, :c))

sys = ssrand(1,1,2)
sys = [sys; sys] # uncontrollable
sysr = minreal(sys)
@test sysr.nx == 2
@test polecompare(sys, [sysr sysr])
@test almost_diagonal(gram(sysr, :c))

# MIMO

sys = ssrand(2,3,4)
sysr = minreal(sys)
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

sys = [sys sys] # unobservable
sysr = minreal(sys)
@test sysr.nx == 4
@test polecompare(sys, [sysr sysr])
@test almost_diagonal(gram(sysr, :c))

sys = ssrand(2,3,4)
sys = [sys; sys] # uncontrollable
sysr = minreal(sys)
@test sysr.nx == 4
@test polecompare(sys, [sysr sysr])
@test almost_diagonal(gram(sysr, :c))


sys = [sys 2sys] 
sysr = minreal(sys)
@test sysr.nx == 4
@test almost_diagonal(gram(sysr, :c))



@testset "diagonalizing" begin
    @info "Testing diagonalizing"
    import ControlSystems: diagonalizing
    P = randn(3,3)
    P = P'P
    T = diagonalizing(P, true)
    @test T*P*T' ≈ I

    T = diagonalizing(P, false)
    TP = T*P*T'
    @test sum(abs, TP - diagm(diag(TP))) < 1e-10

    P = randn(2, 3)
    P = P'P
    T = diagonalizing(P, true)
    @test T*P*T' ≈ ControlSystems.blockdiag(1.0I(2), 0.0I(1))

    T = diagonalizing(P, false)
    TP = T*P*T'
    @test sum(abs, TP - diagm(diag(TP))) < 1e-10

    display(T*P*T')
end

end