using ControlSystems: _equation_order, HammersteinWienerSystem

@testset "equation_order" begin
    @info "Testing equation_order"

    D22 = [0 1; 0 0]
    @test _equation_order(D22) == [2,1]

    D22 = [0 0; 1 0]
    @test _equation_order(D22) == [1,2]

    D22 = [1 0; 0 0]
    @test_throws ErrorException _equation_order(D22)

    function rand_lt(n)
        M = zeros(n, n)
        M[diagind(M, -1)] .= sprandn(n - 1, 0.5)
        M[:, randperm(n)]
        M[diagind(M)] .= 0
        M
    end
    n = 50

    for i = 1:1000
        D = rand_lt(n)
        perm = _equation_order(D)
        @test sort(perm) == 1:n
    end

end



# Test step
println("Simulating first nonlinear system:")
@test step(nonlinearity(identity) * tf(1, [1.0, 1]), 5).y ≈ step(tf(1, [1.0, 1]), 5).y rtol =
    1e-3
@test step(nonlinearity(abs2) * tf(1, [1.0, 1]), 0:0.01:5).y ≈
      abs2.(step(tf(1, [1.0, 1]), 0:0.01:5).y) rtol = 1e-3

# plot(step(nonlinearity(abs2)*tf(1,[1.,1]), 0:0.01:5).y')
# plot!(abs2.(step(tf(1,[1.,1]), 0:0.01:5).y'))

@test step(tf(1, [1.0, 1]) * nonlinearity(x -> 0.5x), 5).y ≈ 0.5step(tf(1, [1.0, 1]), 5).y rtol =
    1e-3

@time step(nonlinearity(abs2) * tf(1, [1.0, 1]))
@time step(nonlinearity(abs2) * tf(1, [1, 1]))
@test step(nonlinearity(abs2) * tf(1, [1, 1])) isa ControlSystemsBase.SimResult


## Multiple nonlinearitys
G = ss(1.0) + nonlinearity(abs2) + nonlinearity(abs2)

s11 = feedback(ss(1 / s), nonlinearity(abs2))
s12 = ss(1 / s)
s21 = HammersteinWienerSystem(ss(1 / (s + 1)))
s22 = ss(1 / (s + 10))

s1 = [s11 s12]
s2 = [s21 s22]

f1 = [s1; s2]
f2 = [
    s11 s12
    s21 s22
]


res = step([s11; s12], 10)
@test size(res.u, 1) == 1
@time y1, t1, x1 = res
@time @test y1[2:2, :] ≈ step(s12, t1)[1] rtol = 1e-14

res = step([s11 s12], 10)
@test size(res.u, 1) == 2
@test size(res.u, 3) == 2

t = 0.0:0.1:10

y1, t1, x1 = step(HammersteinWienerSystem([1.0/s 2/s; 3/s 4/s]), t)
y2, t2, x2 = step([1.0/s 2/s; 3/s 4/s], t)
@test y1 ≈ y2 rtol = 1e-14
@test size(x1, 2) == length(t)
@test size(x1, 3) == 2


## Test only nonlinearity
t = 0:0.01:4
sys = nonlinearity(abs2)

y, t, x = step(sys, t)
@test y[:] ≈ ones(length(t)) atol = 1e-12
@test size(x) == (0, 401, 1)


## Test that nonlinearities can be composed 
sys = nonlinearity(abs2) * nonlinearity(abs2) * 1 / s
res = step(sys, t)

sys2 = nonlinearity(abs2 ∘ abs2) * 1 / s
res2 = step(sys2, t)

@test res.y ≈ res2.y rtol = 1e-3

# plot([res, res2])

## Test predefined nonlinearities
# saturation
using ControlSystemsBase: saturation
th = 0.7
G = tf(1.0, [1, 1])
nl = saturation(th)
C = tf(1)
Cnl = nl * C
L = G * C
Lnl = G * Cnl

@test all(step(feedback(Cnl, G)).y .<= th)

# plot(step([feedback(L); feedback(C, G)], 5), lab = ["Linear y" "Linear u"])
# plot!(step([feedback(Lnl); feedback(Cnl, G)], 5), lab = ["Nonlinear y" "Nonlinear u"])



# offset
using ControlSystemsBase: offset
o = 1.5
G = tf(1.0, [1, 1])
@show nl = offset(o)
@test all(step(nl, 1).y .== 2.5)
# end

w = exp10.(LinRange(-2, 2, 2))
@test all(freqresp(nl, w) .== 1)
@test evalfr(nl, rand())[] == 1

# MIMO offset
o = randn(2)
nl = offset(o)
@test all(step(nl, 1).y[:, :, 1] .== o .+ [1; 0])
@test all(step(nl, 1).y[:, :, 2] .== o .+ [0; 1])

@test all(impulse(nl, 1).y[:, :, 1] .== o)
@test all(impulse(nl, 1).y[:, :, 2] .== o)

G = ssrand(2, 1, 1)
nlG = nl * G
@test nlG.ny1 == G.ny


## Test error on algebraic loop
nl = saturation(0.5)
G = ss(1)

# Two identical saturations in a row is equivalent to a single saturation
# Due to the direct feedthrough in G, we fail to solve this problem and encounter an algebraic loop
D22 = feedback(nl, G*nl).D22
@test_throws ErrorException ControlSystems._equation_order(D22)

## With a proper system, two identical saturations in a row is equivalent to a single saturation
G = ss(-1,1,1,0)
feedback(nl, G*nl)
@test step(feedback(nl, G*nl)).y ≈ step(feedback(nl, G)).y atol=1e-3

## Benchmark
# o = 1.5
# G = tf(1.0, [1, 1])
# nl = offset(o)

# Gnl = feedback(nl*nl*nl*nl*nl*G)
# u = @inline (u,x,t)-> copyto!(u, 0)
# @btime lsim($Gnl, u, 0:0.01:5)
# # 386.459 μs (7695 allocations: 201.34 KiB) @inbounds
# # 395.817 μs (7698 allocations: 201.61 KiB) @inbounds
# # 383.694 μs (7698 allocations: 201.61 KiB) @simd
# # 151.506 μs (1020 allocations: 125.41 KiB) make f a tuple of functions
# # 148.230 μs (1020 allocations: 125.41 KiB) inline u



## Test non-commutativity of multiplication

using ControlSystemsBase: nonlinearity
# build the duffing oscillator 
# ẍ = -kx -k₃x^3 - cẋ + u
# and multiply by 10 on input or output. These are not the same for nonlinear systems
k=10
k3=2
c=1

s = tf("s")
cube = nonlinearity(x->x^3)
vel_loop = feedback(1/s, c)
pos_loop = (k3*cube + k)
duffing = feedback(vel_loop/s, pos_loop)

# Solve for steady state 
# root of (10x + 2x^3) = 10
# https://www.wolframalpha.com/input?i=10x+%2B+2x%5E3+%3D+10
@test step(duffing*10, 1000).y[end] ≈ 0.868830020341475 rtol=1e-3

# 10 times root of (10*x + 2*x^3) = 1
# https://www.wolframalpha.com/input?i=%2810*x+%2B+2*%28x%29%5E3%29+%3D+1
@test step(10*duffing, 1000).y[end] ≈ 10*0.0998011904871354 rtol=1e-3


using ControlSystemsBase
using ControlSystemsBase: ratelimit

l = 0.5
T = 2
for l in 0.1:0.1:1, T in 1:5
    nl = ratelimit(l; Tf=0.01)
    res = step(nl, T)
    # plot(res, plotx=true, plotu=false)
    @test res.y[end] ≈ l*T rtol=0.01
end

##

using ControlSystemsBase
using ControlSystemsBase: deadzone

th = 0.7
@show nl = deadzone(th)
res = lsim(nl, (x,t)->sin(t), 0:0.01:4pi)
# plot(res)
@test minimum(res.y) ≈ -0.3 rtol=1e-3
@test maximum(res.y) ≈ 0.3 rtol=1e-3

@show nl = deadzone(-2th, th)
res = lsim(nl, (x,t)->2sin(t), 0:0.01:4pi)
# plot(res)
@test minimum(res.y) ≈ -0.6 rtol=1e-3
@test maximum(res.y) ≈ 1.3 rtol=1e-3

