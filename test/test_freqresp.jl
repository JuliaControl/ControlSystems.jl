@testset "test_freqresp" begin
## EVALFR ##
H = [tf(0) tf([3, 0],[1, 1, 10]) ; tf([1, 1],[1, 5]) tf([2],[1, 6])]
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])



@test evalfr(H, -6) == [0.0 -0.45; 5.0 Inf]
@test evalfr(H, -5) == [0.0 -0.5; Inf 2.0]
@test evalfr(H, -1) == [0.0 -0.3; 0.0 0.4]
@test evalfr(H, 0) ≈ [0.0 0.0; 0.2 1/3]

@test evalfr(G, -6) == [Inf Inf; Inf Inf]
@test evalfr(G, -5) == [Inf Inf; Inf Inf]
@test evalfr(G, -1) ≈ [0.0 -0.3; 0.0 0.4]
@test evalfr(G, 0) ≈ [0.0 0.0; 0.2 1/3]

## Constant system
w = exp10.(range(-1, stop=1, length=50))

sys1 = ss(1)
G1 = tf(1)
H1 = zpk(1)
resp1 = ones(ComplexF64, length(w), 1, 1)

@test evalfr(sys1, im*w[1]) == fill(resp1[1], 1, 1)
@test evalfr(G1, im*w[1]) == fill(resp1[1], 1, 1)
@test evalfr(H1, im*w[1]) == fill(resp1[1], 1, 1)

@test freqresp(sys1, w) == resp1
@test freqresp(G1, w) == resp1
@test freqresp(H1, w) == resp1

## First order system

sys2 = ss(-1, 1, 1, 1)
G2 = tf([1, 2], [1,1])
H2 = zpk([-2], [-1.0], 1.0)
resp2 = reshape((im*w .+ 2)./(im*w  .+ 1), length(w), 1, 1)

@test evalfr(sys2, im*w[1]) ≈ fill(resp2[1], 1, 1)
@test evalfr(G2, im*w[1]) == fill(resp2[1], 1, 1)
@test evalfr(H2, im*w[1]) == fill(resp2[1], 1, 1)

@test freqresp(sys2, w) ≈ resp2 rtol=1e-15
@test freqresp(G2, w) == resp2
@test freqresp(H2, w) == resp2


## Complex-coefficient system
sys3 = ss(-1+im, 1, (1-im), (1-im))
G3 = (1-im)*tf([1,2-im], [1,1-im])
H3 = zpk([-2+im], [-1+im], (1-im))
resp3 = reshape((1 .- im)*(2 .+ im*(w .- 1))./(1 .+ im*(w .- 1)), length(w), 1, 1)

@test evalfr(sys3, im*w[1]) ≈ fill(resp3[1], 1, 1) rtol=1e-16
@test evalfr(G3, im*w[1]) == fill(resp3[1], 1, 1)
@test evalfr(H3, im*w[1]) == fill(resp3[1], 1, 1)

@test freqresp(sys3, w) ≈ resp3 rtol=1e-15
@test freqresp(G3, w) ≈ resp3 rtol=1e-16
@test freqresp(H3, w) == resp3




## Shortcut notation for evalfr ##
F = tf([1],[1,0.5],-1)
omega = 2
z = 0.5(1+im)
# This test is not correct if we dont have a sample time
# @test_throws F(omega,true)[1] == 1/(exp(-im*2)+0.5)
@test F(z,false)[1] == 1/(z+0.5)
@test_throws ErrorException F(z,true)

F = [tf([1],[1,0.5],2.0) 3*tf([1],[1,0.5],2.0)]
omegas = [1,2]
z = 0.5(1+im)
@test F(omegas[1],true) ≈ [1 3].*1/(exp(im*2)+0.5) atol=1e-14
@test F(omegas[2],true) == [1 3].*1/(exp(2*im*2)+0.5)
@test F(omegas,true) ≈ [k/(exp(omega*im*2)+0.5) for omega=omegas, o=1:1, k=[1,3]] atol=1e-14
@test F(omegas,false) ≈ [k/(omega+0.5) for omega=omegas, o=1:1, k=[1,3]] atol=1e-14
@test F(z,false)[1] == 1/(z+0.5)
@test_throws ErrorException F(z,true)

## Test bode, nyquist and sigma
sys = [tf([1,-1], [1,1,1]) 0; 0 tf([1],[1,1])]
f(s) = [(s-1)./(s.^2+s+1) 0; 0 1 ./(1+s)]
ws = exp10.(range(-2, stop=2, length=50))
resp = Array{ComplexF64}(undef, 50,2,2)
for (i,w) in enumerate(ws)
    resp[i,:,:] = f(im*w)
end

@test bode(sys, ws)[1:2] == (abs.(resp), rad2deg.(angle.(resp)))
@test nyquist(sys, ws)[1:2] == (real(resp), imag(resp))
sigs = Array{Float64}(undef, 50,2)
for i in eachindex(ws)
    sigs[i,:] =  svdvals(resp[i,:,:])
end
@test sigma(sys, ws)[1] == sigs

#Test default freq vector contains at least half a decade more than all features
p, z = 100, 1/1000
sys2 = [tf([1],[1/p,1]) tf([1/z, 2/z, 1],[1])]
mag, phase, ws2 = bode(sys2)
@test maximum(ws2) >= 2max(p,z)
@test minimum(ws2) <= 0.5min(p,z)

#Test default freq vector not too small
p, z = 1, 1.2
sys2 = [tf([1],[1/p,1]) tf([1/z, 2/z, 1],[1])]
mag, mag, ws2 = bode(sys2)
@test maximum(ws2) >= 5max(p,z)
@test minimum(ws2) <= 0.2min(p,z)
@test length(ws2) > 100
end
