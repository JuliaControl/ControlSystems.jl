import DelayDiffEq: MethodOfSteps, Tsit5
using SparseArrays

s = tf('s')
s11 = feedback(ss(1/s), delay(1))
s12 = ss(1/s)
s21 = DelayLtiSystem(ss(1/(s+1)))
s22 = ss(1/(s+10))

s1 = [s11 s12]
s2 = [s21 s22]

f1 = [s1;s2]
f2 = [s11 s12;
      s21 s22]


# Test step
println("Simulating first delay system:")
@time step(delay(1)*tf(1,[1.,1]))
@time step(delay(1)*tf(1,[1,1]))
@test step(delay(1)*tf(1,[1,1])) isa ControlSystemsBase.SimResult

res = step([s11;s12], 10)
@test size(res.u,1) == 1
@time y1, t1, x1 = res
@time @test y1[2:2,:] ≈ step(s12, t1)[1] rtol = 1e-14

res = step([s11 s12], 10)
@test size(res.u,1) == 2
@test size(res.u,3) == 2

t = 0.0:0.1:10
y2, t2, x2 = step(s1, t)
# TODO Figure out which is inexact here
@test y2[1:1,:,1:1] + y2[1:1,:,2:2] ≈ step(s11, t)[1] + step(s12, t)[1] rtol=1e-14

y3, t3, x3 = step([s11; s12], t)
@test y3[1:1,:,1] ≈ step(s11, t)[1] rtol = 1e-14
@test y3[2:2,:,1] ≈ step(s12, t)[1] rtol = 1e-14

y1, t1, x1 = step(DelayLtiSystem([1.0/s 2/s; 3/s 4/s]), t)
y2, t2, x2 = step([1.0/s 2/s; 3/s 4/s], t)
@test y1 ≈ y2 rtol=1e-14
@test size(x1,2) == length(t)
@test size(x1,3) == 2



##################### Test with known solution
K = 2.0
sys_known = feedback(delay(1)*K/s, 1)

ystep, t, _ = step(sys_known, 3)

function y_expected(t, K)
      if t < 1
            return 0
      elseif t < 2
            return K*(t-1)
      elseif t <= 3
            return K*(t-1)-1/2*K^2*(t-2)^2
      else
            throw(ArgumentError("Test not defined here"))
      end
end

@test ystep' ≈ y_expected.(t, K) atol = 1e-12

function dy_expected(t, K)
      if t < 1
            return 0
      elseif t < 2
            return K
      elseif t <= 3
            return K - K^2*(t-2)
      else
            throw(ArgumentError("Test not defined here"))
      end
end

y_impulse, t, _ = impulse(sys_known, 3)

# TODO Better accuracy for impulse
@test y_impulse' ≈ dy_expected.(t, K) rtol=1e-13
@test maximum(abs, y_impulse' - dy_expected.(t, K)) < 1e-12

y_impulse, t, _ = impulse([sys_known sys_known], 3)

@test y_impulse[1,:,1] ≈ dy_expected.(t, K) rtol=1e-13
@test y_impulse[1,:,2] ≈ dy_expected.(t, K) rtol=1e-13

y_impulse, t, _ = impulse(sys_known, 3, alg=MethodOfSteps(Tsit5()))
# Two orders of magnitude better with BS3 in this case, which is default for impulse
@test y_impulse' ≈ dy_expected.(t, K) rtol=1e-5
@test maximum(abs, y_impulse' - dy_expected.(t, K)) < 1e-5

## Test delay with D22 term
t = 0:0.01:4

sys = delay(1)

y, t, x = step(sys, t)
@test y[:] ≈ [zeros(100); ones(301)] atol = 1e-12
@test size(x) == (0,401)

sys = delay(1)*delay(1)*1/s

y, t, x = step(sys, t)

y_sol = [zeros(200);0:0.01:2]'

@test maximum(abs,y-y_sol) < 1e-13
@test maximum(abs,x-collect(0:0.01:4)') < 1e-15

# TODO For some reason really bad accuracy here
# Looks like a lag in time
sys = 1/s*delay(1)*delay(1)

y, t, x = step(sys, t)
@test maximum(abs,y-y_sol) < 1e-5
@test maximum(abs,x-y_sol) < 1e-5

t = 0:0.001:0.1
y, t, x = step(sys, t)
@test length(y) == length(t)


# Step below fails without force_dtmin to the solver
s = tf("s")
P = 1 / (0.85*s + 1)*exp(-0.14*s)
res = step(P, 5)
@test res.t[end] > 4.5
@test length(res.y) > 30