# Tuning a PID controller from data

In this example, we will consider a very commonly occurring workflow: using process data to tune a PID controller.

The two main steps involved in this workflow are:
1. Estimate a process model from data
2. Design a controller based on the estimated model

In this example, which is split into two parts, we will consider tuning a velocity controller for a **flexible robot arm**. Part 1 is available here: [Flexible Robot Arm Part 1: Estimation of a model.](https://baggepinnen.github.io/ControlSystemIdentification.jl/dev/examples/flexible_robot/). The system identification uses the package [ControlSystemIdentification.jl](https://github.com/baggepinnen/ControlSystemIdentification.jl).

The rest of this example makes up part 2, tuning of the controller. We simply replicate the relevant code from part 1 to get the estimated model, and then use the estimated model to tune controllers.
```@example PID_TUNING
using DelimitedFiles, Plots
using ControlSystemIdentification, ControlSystems

url = "https://ftp.esat.kuleuven.be/pub/SISTA/data/mechanical/robot_arm.dat.gz"
zipfilename = "/tmp/flex.dat.gz"
path = Base.download(url, zipfilename)
run(`gunzip -f $path`)
data = readdlm(path[1:end-3])
u = data[:, 1]' # torque
y = data[:, 2]' # acceleration
d = iddata(y, u, 0.01) # sample time not specified for data, 0.01 is a guess
Pacc = subspaceid(d, 4, focus=:prediction) # Estimate the process model using subspace-based identification
```


Since the data used for the system identification had acceleration rather than velocity as output, we multiply the estimated model by the transfer function ``1/s`` to get a velocity model. Before we do this, we convert the estimated discrete-time model into continuous time using the function [`d2c`](@ref). The estimated system also has a negative gain due to the mounting of the accelerometer, so we multiply the model by ``-1`` to get a positive gain.
```@example PID_TUNING
s = tf("s")
P = 1/s * d2c(-Pacc.sys)
bodeplot(P)
```

## Controller tuning
We could take multiple different approaches to tuning the PID controller, a few alternatives are listed here
- Trial and error in simulation or experiment.
- Manual loop shaping
- Automatic loop shaping
- Step-response optimization ([example](https://juliacontrol.github.io/ControlSystems.jl/stable/examples/automatic_differentiation/#Optimization-based-tuning%E2%80%93PID-controller))

Here, we will attempt a manual loop-shaping approach using the function [`loopshapingPID`](@ref), and then then compare the result to a pole-placement controller.

### Manual loop shaping
The function [`loopshapingPID`](@ref) takes a model and selects the parameters of a PID-controller such that the Nyquist curve of the loop-transfer function ``L = PC`` at the frequency `ω` is tangent to the circle where the magnitude of the complimentary sensitivity function ``T = PC / (1+PC)`` equals ``M_T``. This allows us to explicitly solve for the PID parameters that achieves a desired target value of ``M_T`` at a desired frequency `ω`. The function can optionally produce a plot which draws the design characteristics and the resulting Nyquist curve, which we will make use of here. A Youtube video tutorial that goes into more details on how this function works is [available here](https://youtu.be/BolNmqYYIEg?si=hF-xvsPL_wBngpft&t=775).

Since the process contains two sharp resonance peaks, visible in the Bode diagram above, the requirements for our velocity controller have to be rather modest. We therefore tell [`loopshapingPID`](@ref) that we want to include a lowpass filter in the controller to suppress any frequencies above ``ω_f = 1/T_f`` so that the resonances do not cause excessive robustness problems. We choose the design frequency to be ``ω = 5`` and the target value of ``M_T = 1.35`` achieved at an angle of ``ϕ_t = 35`` degrees from the negative real axis. The function returns the controller, the PID parameters, the resulting Nyquist curve, and the lowpass-filter controlled `CF`. 

```@example PID_TUNING
ω = 5
Tf = 1/10
C, kp, ki, kd, fig, CF = loopshapingPID(P, ω; Mt = 1.35, ϕt=35, doplot=true, Tf)
fig
```
The PID parameters are by default returned on "standard form", but the parameter convention to use can be selected using the `form` keyword.

The result above satisfies the design in the design point, but the sharp resonances peak well above the desired maximum of the complementary sensitivity function. The problem here is that a PID controller is fundamentally incapable at damping the resonances in this high-order system. Indeed, we have a closed-loop system with a 8-dimensional state, but only 3-4 parameters in the PID controller (depending on whether or not we count the filter parameter), so there is no hope for us to arbitrarily place the poles using the PID controller. This can result in poor robustness properties, as we will see below.

Next, we form the closed-loop system ``G`` from reference to output an plot a step response
```@example PID_TUNING
G = feedback(P*CF)
plot(step(G, 10), label="Step response")
```
This looks extremely aggressive and with clear resonances visible. The problem here is that no mechanical system can follow a perfect step in the reference, and it is thus common to generate some form of physically realizable smooth step as input reference. Below, we use the package [TrajectoryLimiters.jl](https://github.com/baggepinnen/TrajectoryLimiters.jl) to filter the reference step such that it has bounded acceleration and velocity

```@example PID_TUNING
using TrajectoryLimiters
ẋM = 2 # Velocity limit
ẍM = 1 # Acceleration limit
limiter = TrajectoryLimiter(d.Ts, ẋM, ẍM)
inputstep, vel, acc = limiter([0; ones(1000)])
timevec = 0:d.Ts:10
plot(step(G, 10), label="Step response")
plot!(lsim(G, inputstep', timevec), label="Smooth step response")
plot!(timevec, inputstep, label="Smooth reference trajectory", l=(:dash, :black))
```

The result now looks much better, with some small amount of overshoot. The performance is not terrific, taking about 2 seconds to realize the step. However, attempting to make the response faster using feedback alone will further exacerbate the robustness problems due to the resonance peaks highlighted above.

A more conservative and robust tuning, that does not let the resonance peaks cause large peaks in the sensitivity functions, can be realized by manual loop shaping with the help of a [`marginplot`](@ref)
```@example PID_TUNING
Tf = 0.4
Ti = 4
Td = 0.1
CF = pid(10, Ti, Td; Tf)
marginplot(P*CF)
```

This tuning shows good gain and phase margins, but the price we pay for this is of course performance:
```@example PID_TUNING
ẍM = 0.008 # Acceleration limit
limiter2 = TrajectoryLimiter(d.Ts, ẋM, ẍM)
inputstep2, vel, acc = limiter2([0; ones(5000)])
timevec = 0:d.Ts:50
G = feedback(P*CF)
plot(step(G, 50), label="Step response")
plot!(lsim(G, inputstep2', timevec), label="Smooth step response")
plot!(timevec, inputstep2, label="Smooth reference trajectory", l=(:dash, :black))
```
The closed-loop system now responds significantly slower. 

Below, we attempt a pole-placement design for comparison. Contrary to the PID controller, a pole-placement controller _can_ place all poles of this system arbitrarily (the system is _controllable_, which can be verified using the function [`controllability`](@ref)).


## Pole placement
We start by inspecting the pole locations of the open-loop plant
```@example PID_TUNING
pzmap(P)
```
As expected, we have 2 resonant pole pairs.

When dampening fast resonant poles, it is often a good idea to _only_ dampen them, not to change the bandwidth of them. Trying to increase the bandwidth of these fast poles requires very large controller gain, and making the poles slower often causes severe robustness problems. We thus place the resonant poles with the same magnitude, but with perfect damping.
```@example PID_TUNING
current_pole_magnitudes = abs.(poles(P))
```

The integrator pole can be placed to achieve a desired bandwidth. Here, we place it in -25rad/s to achieve a faster response than the PID controller achieved.
```@example PID_TUNING
desired_poles = -[80, 80, 37, 37, 25]
```

We compute the state-feedback gain ``L`` using the function [`place`](@ref), and also compute an observer gain ``K`` using the rule of thumb that the observer poles should be approximately twice as fast as the system poles.
```@example PID_TUNING
L = place(P, desired_poles, :c)
K = place(P, 2*desired_poles, :o)
```

The resulting observer-based state-feedback controller can be constructed using the function [`observer_controller`](@ref). We also form the closed-loop system ``G_{pp}`` from reference to output an plot a step response like we did above
```@example PID_TUNING
Cpp = observer_controller(P, L, K)
Gpp = feedback(P*Cpp)
plot(lsim(Gpp, inputstep', timevec), label="Smooth step response")
plot!(timevec, inputstep, label="Smooth reference trajectory")
```
The pole-placement controller achieves a very nice result, but this comes at a cost of using very large controller gain. The gang-of-four plot below indicates that we have a controller with reasonable robustness properties if we inspect the sensitivity and complimentary sensitivity functions, but the noise-amplification transfer function ``CS`` has a large gain for high frequencies, implying that this controller requires a very good sensor to be practical!
```@example PID_TUNING
gangoffourplot(P, Cpp)
```

With the PID controller, we can transform the PID parameters to the desired form and enter those into an already existing PID-controller implementation. Care must be taken to incorporate also the measurement filter designed by [`loopshapingPID`](@ref), this filter is important for robustness analysis to be valid. If no existing PID controller implementation is available, we may either make use of the package [DiscretePIDs.jl](https://github.com/JuliaControl/DiscretePIDs.jl), or generate C-code for the controller. Below, we generate some C code.


## C-Code generation
Using the pole-placement controller derived above, we discretize the controller using the Tustin (bilinear) method of the function [`c2d`](@ref), and then call [`SymbolicControlSystems.ccode`](https://github.com/JuliaControl/SymbolicControlSystems.jl#code-generation).
```julia
using SymbolicControlSystems
Cdiscrete = c2d(Cpp, d.Ts, :tustin)
SymbolicControlSystems.ccode(Cdiscrete)
```

This produces the following C-code for filtering the error signal through the controller transfer function
```c
#include <stdio.h>

#include <math.h>

void transfer_function(double *y, double u) {
    static double x[5] = {0};  // Current state

    double xp[5] = {0};        // Next state
    int i;

    // Advance the state xp = Ax + Bu
    xp[0] = (1.323555302697655*u - 0.39039743126198218*x[0] - 0.0016921457205018749*x[1] - 0.0012917116898466163*x[2] + 0.001714187010327197*x[3] + 0.0016847122113737578*x[4]);
    xp[1] = (96.429820608571958*u - 95.054670090613683*x[0] + 0.13062589956122247*x[1] + 0.78522537641468981*x[2] - 0.21646419099004577*x[3] - 0.081049292550184435*x[4]);
    xp[2] = (13.742733359914924*u - 15.008953114410946*x[0] - 0.89468526010523608*x[1] + 0.717920592086567*x[2] + 0.025588437849127441*x[3] + 0.021322717715438085*x[4]);
    xp[3] = (303.50179259195619*u - 268.71085904944562*x[0] + 1.251906632298234*x[1] + 0.62490471615521814*x[2] + 0.15988074336172073*x[3] - 0.4891888301891486*x[4]);
    xp[4] = (-27.542490469297601*u + 37.631007484177218*x[0] + 1.2366332766644277*x[1] + 0.11855488877285068*x[2] - 0.29543727245267387*x[3] + 0.76660106988104448*x[4]);

    // Accumulate the output y = C*x + D*u
    y[0] = (912.01950044640216*u - 742.76679702406477*x[0] + 10.364451210258789*x[1] + 2.7824392013821053*x[2] - 2.3907024395896719*x[3] - 3.734615363051947*x[4]);

    // Make the predicted state the current state
    for (i=0; i < 5; ++i) {
        x[i] = xp[i];
    }

}
```

## Summary
This tutorial has shown how to follow a workflow that consists of
1. Estimate a process model using experimental data.
2. Design a controller based on the estimated model.
3. Simulate the closed-loop system and analyze its robustness properties.
4. Generate C-code for the controller.

Each of these steps is covered in additional detail in the videos available in the playlist [Control systems in Julia](https://youtube.com/playlist?list=PLC0QOsNQS8hZtOQPHdtul3kpQwMOBL8Qc&si=yUrXz5cH4QqTPlR_). See also the tutorial [Control design for a quadruple-tank system](https://help.juliahub.com/juliasimcontrol/dev/examples/quadtank/).