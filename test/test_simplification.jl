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
@test minreal(sys-sysmin, atol=eps()) == ss(0)

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
@test tf(sys) ≈ tf(sysr) # no change
@test hinfnorm(sys - sysr)[1] < sqrt(eps())
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

syse = [sys sys] # unobservable
sysr = minreal(syse)
@test tf(sys) ≈ tf(sysr)
@test hinfnorm(syse - sysr)[1] < sqrt(eps())
@test sysr.nx == 2
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

sys = ssrand(1,1,2)
syse = [sys; sys] # uncontrollable
sysr = minreal(syse)
@test tf(sys) ≈ tf(sysr)
@test hinfnorm(syse - sysr)[1] < sqrt(eps())
@test sysr.nx == 2
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

# MIMO

sys = ssrand(2,3,4)
sysr = minreal(sys)
@test tf(sys) ≈ tf(sysr)
@test hinfnorm(sys - sysr)[1] < sqrt(eps())
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

syse = [sys sys] # unobservable
sysr = minreal(syse)
@test tf(sys) ≈ tf(sysr)
@test hinfnorm(syse - sysr)[1] < sqrt(eps())
@test sysr.nx == 4
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))

sys = ssrand(2,3,4)
syse = [sys; sys] # uncontrollable
sysr = minreal(syse)
@test hinfnorm(syse - sysr)[1] < sqrt(eps())
@test sysr.nx == 4
@test polecompare(sys, sysr)
@test almost_diagonal(gram(sysr, :c))


syse = [syse 2syse] 
sysr = minreal(syse)
@test hinfnorm(syse - sysr)[1] < sqrt(eps())
@test sysr.nx == 4
@test almost_diagonal(gram(sysr, :c))

## A difficult test case
A = [0 0.0078125 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0.03125 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0.0625 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0.125 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0.25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0.25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0.25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0.25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4.0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8.0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 16; -0.011884937660980725 -0.012778154166964564 -0.01119260178557217 -0.009637610244725788 -0.013444333518321703 -0.024933403383717623 -0.08555297518983913 -0.10547615432200363 -0.257315496122918 -0.23187549363330312 -0.21667524574733332 -0.30223256771695023 -0.45121173605407 -0.5047893604121784 -0.6178244657434898 -0.5670334792048756 -0.578825622781816 -0.44204200198750454 -0.761671019791637 -0.9760854687720177 -1.43125569071897 -1.544690913834306 -1.938382943806063 -1.761496777281267 -1.898323520495274 -1.446578849898936 -1.3410048784833273 -0.8495997571966917 -1.3547082946323532 -1.4069711503743814 -1.9251377271914922 -1.602558738496262 -1.8726415465003627 -1.2061444209280585 -1.1935187030477974 -1.1218409477943014 -1.853591287718073 -1.1389270252587735 -3.0706102237822135 -3.751580204015136]
B = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2.0]
C = [0.005812661478252773 0.006433233648168637 0.005569209729056593 0.0048397954489543486 0.006663369378465938 0.012431321785229477 0.04193362776128268 0.051873916417509304 0.12367234268488136 0.11157929309057098 0.1011015880191636 0.14095709673373336 0.20204618897580154 0.22557184889576087 0.2618906601439306 0.2393828503415872 0.22851174487610879 0.17328083467568964 0.2746761928818859 0.3478917988051352 0.460864762256988 0.4881063153744335 0.5428129947045572 0.47905729320817864 0.448626601273217 0.3272125013158462 0.25850394937566856 0.1537263533721266 0.20483555313250637 0.19464007594416033 0.21782566127726125 0.1604008091677848 0.14918002448754064 0.08114818712985343 0.06131756636534406 0.04536653837340357 0.053073352879493864 0.022407590502674234 0.03575549467417118 0.017405094415233746]
D = [0.0]
sys = ss(A,B,C,D)
sysr = minreal(sys)
@test hinfnorm(sys - sysr)[1] < sqrt(eps())



@testset "diagonalizing" begin
    @info "Testing diagonalizing"
    import ControlSystems: hermitian_diagonalizing
    P = randn(3,3)
    P = P'P
    T = hermitian_diagonalizing(P, true)
    @test T*P*T' ≈ I

    T = hermitian_diagonalizing(P, false)
    TP = T*P*T'
    @test sum(abs, TP - diagm(diag(TP))) < 1e-10

    P = randn(2, 3)
    P = P'P # create rank-deficient matrix of rank 2
    T = hermitian_diagonalizing(P, true)
    @test T*P*T' ≈ ControlSystems.blockdiag(1.0I(2), 0.0I(1)) # should result in one 0 on diag

    T = hermitian_diagonalizing(P, false)
    TP = T*P*T'
    @test sum(abs, TP - diagm(diag(TP))) < 1e-10

    display(T*P*T')
end

end