@testset "test_freqresp" begin
## EVALFR ##
H = [tf(0) tf([3, 0],[1, 1, 10]) ; tf([1, 1],[1, 5]) tf([2],[1, 6])]
G = ss([-5 0 0 0; 0 -1 -2.5 0; 0 4 0 0; 0 0 0 -6], [2 0; 0 1; 0 0; 0 2],
       [0 3 0 0; -2 0 0 1], [0 0; 1 0])



@test evalfr(H, -6) == [0.0 -0.45; 5.0 Inf]
@test evalfr(H, -5) == [0.0 -0.5; Inf 2.0]
@test evalfr(H, -1) == [0.0 -0.3; 0.0 0.4]
@test evalfr(H, 0) ≈ [0.0 0.0; 0.2 1/3]
@inferred evalfr(H, 0)

@test (@test_logs (:warn, "Got exception SingularException(4), returning Inf") evalfr(G, -6)) == [Inf Inf; Inf Inf]
@test (@test_logs (:warn, "Got exception SingularException(1), returning Inf") evalfr(G, -5)) == [Inf Inf; Inf Inf]
@test evalfr(G, -1) ≈ [0.0 -0.3; 0.0 0.4]
@test evalfr(G, 0) ≈ [0.0 0.0; 0.2 1/3]

## Constant system
w = exp10.(range(-1, stop=1, length=50))

sys1 = ss(1)
sys1s = HeteroStateSpace(sparse([1]))
G1 = tf(1)
H1 = zpk(1)
resp1 = ones(ComplexF64, length(w), 1, 1)

@test evalfr(sys1, im*w[1]) == fill(resp1[1], 1, 1)
@test evalfr(G1, im*w[1]) == fill(resp1[1], 1, 1)
@test evalfr(H1, im*w[1]) == fill(resp1[1], 1, 1)

@test freqresp(sys1, w) == resp1
@test freqresp(sys1s, w) == resp1
@test freqresp(G1, w) == resp1
@test freqresp(H1, w) == resp1
@inferred freqresp(sys1, w)
@inferred freqresp(G1, w)
@inferred freqresp(H1, w)

for G in [2, 2I, 2I(2), randn(2,2)]
    @test dcgain(G) == G
    !(G isa UniformScaling) && @test freqresp(G, [1,2]) == cat(G,G, dims=3) # I does not have size
    @test freqresp(G, 1) == G
    @test evalfr(G, 1im) == G
end

## First order system

sys2 = ss(-1, 1, 1, 1)
sys2s = HeteroStateSpace(sparse([-1;;]),sparse([1;;]),sparse([1;;]),sparse([1;;])) # test that freqresp works for matrix types that don't support hessenberg
G2 = tf([1, 2], [1,1])
H2 = zpk([-2], [-1.0], 1.0)
resp2 = reshape((im*w .+ 2)./(im*w  .+ 1), length(w), 1, 1)

@test evalfr(sys2, im*w[1]) ≈ fill(resp2[1], 1, 1)
@test evalfr(G2, im*w[1]) == fill(resp2[1], 1, 1)
@test evalfr(H2, im*w[1]) == fill(resp2[1], 1, 1)

@test freqresp(sys2, w) ≈ resp2 rtol=1e-15
@test freqresp(sys2s, w) ≈ resp2 rtol=1e-15
@test freqresp(G2, w) == resp2
@test freqresp(H2, w) == resp2

@inferred freqresp(sys2, w)
@inferred freqresp(sys2s, w)
@inferred freqresp(G2, w)
@inferred freqresp(H2, w)
@test (@allocated freqresp(sys2, w)) < 1.2*56672 # allow 20% increase due to compiler variations
@test (@allocated freqresp(G2, w)) < 1.2*976 # allow 20% increase due to compiler variations

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


## Benchmark code for freqresp
# sizes = [1:40; 50:5:100; 120:20:300; 800]
# times1 = map(sizes) do nx
#     w = exp10.(LinRange(-2, 2, 200))
#     @show nx
#     G = ssrand(2,2,nx)
#     nx == 1 && (freqresp(G, w)) # precompile
#     GC.gc()
#     t = @timed freqresp(G, w)
#     (t.time, t.bytes)
# end
# sleep(5)
# times2 = map(sizes) do nx
#     w = exp10.(LinRange(-2, 2, 200))
#     @show nx
#     G = ssrand(2,2,nx)
#     # NOTE: rename the freqresp method to be benchmarked before running. Below it's called freqresp_large
#     nx == 1 && (freqresp_large(G, w)) # precompile
#     GC.gc()
#     t = @timed freqresp_large(G, w)
#     (t.time, t.bytes)
# end

# f1 = plot(sizes, first.(times1), scale=:log10, lab="Time freqresp", m=:o)
# plot!(sizes, first.(times2), scale=:log10, lab="Time freqresp_large", xlabel="Model order", m=:o)

# f2 = plot(sizes, last.(times1), scale=:log10, lab="Allocations freqresp", m=:o)
# plot!(sizes, last.(times2), scale=:log10, lab="Allocations freqresp_large", xlabel="Model order", m=:o)
# plot(f1, f2)