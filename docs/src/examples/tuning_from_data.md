# Tuning a PID controller from data

In this example, we will consider a very commonly occurring workflow: using process data to tune a PID controller.

The two main steps involved in this workflow are:
1. Estimate a process model from data
2. Characterize the uncertainty in the estimated model
3. Design a controller based on the estimated model
4. Verify that the controller is robust with respect to the estimated model uncertainty


## Estimation of a model
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



## Dealing with model uncertainty
When using a model for control design, we always have to consider how robust we are with respect to errors in the model. Classical margins like the gain and phase margins are simple measures of robustness, applicable to simple measures of uncertainty. Here, we will attempt to characterize the uncertainty in the model slightly more accurately.

When we estimate linear black-box models from data, like we did above using `subspaceid`, we can get a rough estimate of how well a linear model describes the input-output data by looking at the magnitude-squared coherence function ``\gamma(i\omega)``:
```@example PID_TUNING
coherenceplot(d)
```
For frequencies where ``\gamma`` is close to one, a linear model is expected to fit well, whereas for frequencies where ``\gamma`` is close to zero, we cannot trust the model. How does this rough estimate of model certainty translate to our control analysis? In the video [The benefit and Cost of Feedback](https://youtu.be/uQx192FyA5g?si=kubWnq__ohWOaICw), we show that for frequencies where the uncertainty in the model is large, we must have a small sensitivity. In the video, we analyzed the effects of additive uncertainty, in which case we need to make sure that the sensitivity function ``CS = C/(1+PC)`` is sufficiently small. When using the rough estimate of model uncertainty provided by the coherence function, it may be more reasonable to consider a multiplicative (relative) uncertainty model, in which case we need to verify that the sensitivity function ``T = PC/(1+PC)`` is small for frequencies where ``\gamma`` is small.

Since our coherence drops significantly above ``\omega = 130``rad/s, we will try to design a controller that yields a complementary sensitivity function ``T`` that has low gain above this frequency.


In the [documentation of RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/uncertainty/#Multiplicative-uncertainty), we list a number of common uncertainty models together with the criteria for robust stability. A good resource on the gang of four is available in [these slides](https://www.control.lth.se/fileadmin/control/staff/KJ/FeedbackFundamentals.pdf).

## Controller tuning
We could take multiple different approaches to tuning the PID controller, a few alternatives are listed here
- Trial and error in simulation or experiment.
- Manual loop shaping
- Automatic loop shaping
- Step-response optimization ([example](https://juliacontrol.github.io/ControlSystems.jl/stable/examples/automatic_differentiation/#Optimization-based-tuning%E2%80%93PID-controller))

Here, we will attempt a manual loop-shaping approach using the function [`loopshapingPID`](@ref), and then then compare the result to a pole-placement controller.

### Manual loop shaping


Since the process contains two sharp resonance peaks, visible in the Bode diagram above, we want to include a lowpass filter in the controller to suppress any frequencies above the first resonance so that the resonances do not cause excessive robustness problems. Here, we will use a second-order lowpass filer.

A PID controller is fundamentally incapable at damping the resonances in this high-order system. Indeed, for a plant model of order 4, we have a closed-loop system with a 7-dimensional state (one pole for the integrator and two for the low-pass filter), but only 3-4 parameters in the PID controller (depending on whether or not we count the filter parameter), so there is no hope for us to arbitrarily place the poles using the PID controller. Trying to use a gain high enough to dampen the resonant poles can result in poor robustness properties, as we will see below.


The function [`pid`](@ref) takes the PID parameters "standard form", but the parameter convention to use can be selected using the `form` keyword. We use the function [`marginplot`](@ref) to guide our tuning, the following parameters were found to give a good result
```@example PID_TUNING
K = 10
Tf = 0.4
Ti = 4
Td = 0.1
CF = pid(K, Ti, Td; Tf)
marginplot(P*CF)
```
Here, we have selected the proportional gain ``K``large enough to give a crossover bandwidth of about 1rad/s, being careful not to let the resonance peaks reach too close to unit gain, destroying our robustness. The integral time constant ``T_i`` is selected as low as possible without destroying the phase margin, and the derivative time constant ``T_d`` is increased slowly to improve the phase margin while not letting the resonance peaks become too large.

The [`pid`](@ref) function returns the PI controller with the second-order lowpass filter already applied.

Next, we form the closed-loop system ``G`` from reference to output an plot a step response
```@example PID_TUNING
G = feedback(P*CF)
plot(step(G, 50), label="Step response")
```
This looks rather aggressive and with a large overshoot visible. The problem here is that no mechanical system can follow a perfect step in the reference, and it is thus common to generate some form of physically realizable smooth step as input reference. Below, we use the package [TrajectoryLimiters.jl](https://github.com/baggepinnen/TrajectoryLimiters.jl) to filter the reference step such that it has bounded acceleration and velocity

```@example PID_TUNING
using TrajectoryLimiters
ẋM = 1 # Velocity limit
ẍM = 0.01 # Acceleration limit
limiter = TrajectoryLimiter(d.Ts, ẋM, ẍM)
inputstep, vel, acc = limiter([0; ones(5000)])
timevec = 0:d.Ts:50
plot(step(G, 50), label="Step response")
plot!(lsim(G, inputstep', timevec), label="Smooth step response")
plot!(timevec, inputstep, label="Smooth reference trajectory", l=(:dash, :black))
```

The result now looks much better, with some small amount of overshoot. The performance is not terrific, taking about 20 seconds to realize the step. However, attempting to make the response faster using feedback alone will further exacerbate the robustness problems due to the resonance peaks highlighted above.



To analyze the robustness of this controller, we can inspect the sensitivity functions in the [gang of four](https://www.control.lth.se/fileadmin/control/staff/KJ/FeedbackFundamentals.pdf). In particular, we are interested in the complementary sensitivity function ``T = PC/(1+PC)`` 
```@example PID_TUNING
gangoffourplot(P, CF)
```

The gang of four indicates that we have a robust tuning, no uncomfortably large peaks appears in either ``T`` or ``S``.

Below, we attempt a pole-placement design for comparison. Contrary to the PID controller, a pole-placement controller _can_ place all poles of this system arbitrarily (the system is _controllable_, which can be verified using the function [`controllability`](@ref)).


## Pole placement
We start by inspecting the pole locations of the open-loop plant
```@example PID_TUNING
pzmap(P)
```
As expected, we have 2 resonant pole pairs.

When dampening fast resonant poles, it is often a good idea to _only_ dampen them, not to change the bandwidth of them. Trying to increase the bandwidth of these fast poles requires very large controller gain, and making the poles slower often causes severe robustness problems. We thus try to place the resonant poles with the same magnitude, but with perfect damping.
```@example PID_TUNING
current_poles = poles(P)
```

The integrator pole can be placed to achieve a desired bandwidth. Here, we place it in -30rad/s to achieve a faster response than the PID controller achieved.
```@example PID_TUNING
desired_poles = -[80, 80, 37, 37, 30];
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
plot!(timevec, inputstep, label="Smooth reference trajectory", l=(:dash, :black), legend=:bottomright)
```
The pole-placement controller achieves a very nice result, but this comes at a cost of using very large controller gain. The gang-of-four plot below indicates that we have a controller with a very large noise-amplification transfer function, ``CS`` has a large gain for high frequencies, implying that this controller requires a very good sensor to be practical! We also have significant gain in ``T`` well above the frequency ``ω = 130``rad/s above which we couldn't trust the model.
```@example PID_TUNING
gangoffourplot(P, Cpp)
vline!(fill(130, 1, 4), label="\$ω = 130\$", l=(:dash, :black))
```


Due to the high gain of the controller we got, we redo the design, this time only dampening the resonant poles slightly. We also lower the bandwidth of the integrator pole to make the controller less aggressive
```@example PID_TUNING
p1 = current_poles[2]
p2 = current_poles[4]

p1_new = abs(p1) * cis(-pi + deg2rad(65)) # Place the pole with the same magnitude, but with an angle of -pi + 65 degrees
p2_new = abs(p2) * cis(-pi + deg2rad(65))
desired_poles = [-20, p1_new, conj(p1_new), p2_new, conj(p2_new)]
L = place(P, desired_poles, :c) |> real
K = place(P, 2*desired_poles, :o) |> real
Cpp = observer_controller(P, L, K)
Gpp = feedback(P*Cpp)
f1 = plot(lsim(Gpp, inputstep', timevec), label="Smooth step response")
plot!(timevec, inputstep, label="Smooth reference trajectory", l=(:dash, :black), legend=:bottomright)

f2 = gangoffourplot(P, Cpp)
vline!(fill(130, 1, 4), label="\$ω = 130\$", l=(:dash, :black))
plot(f1, f2, size=(800, 600))
```
We still have a nice step response using this controller, but this time, we have a rolloff in ``T`` that starts around the frequency ``ω = 130``rad/s.

## C-Code generation
With the PID controller, we can transform the PID parameters to the desired form and enter those into an already existing PID-controller implementation. Care must be taken to incorporate also the measurement filter designed, this filter is important for robustness analysis to be valid. If no existing PID controller implementation is available, we may either make use of the package [DiscretePIDs.jl](https://github.com/JuliaControl/DiscretePIDs.jl), or generate C-code for the controller. Below, we generate some C code.


Using the pole-placement controller derived above, we discretize the controller using the Tustin (bilinear) method with the function [`c2d`](@ref), and then call [`SymbolicControlSystems.ccode`](https://github.com/JuliaControl/SymbolicControlSystems.jl#code-generation) to generate the code.
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
    xp[0] = (1.2608412916795442*u - 0.35051915762703334*x[0] + 0.0018847792079810998*x[1] - 0.0035104037080211504*x[2] + 0.0022503125347378308*x[3] + 0.00019318421187795658*x[4]);
    xp[1] = (45.346976964169166*u - 49.856146529754966*x[0] + 0.19058339536496746*x[1] + 0.58214123400704609*x[2] - 0.068048140252114517*x[3] - 0.03667586076286556*x[4]);
    xp[2] = (18.14135831274827*u - 19.16237014106056*x[0] - 0.84117137404200237*x[1] + 0.7024229589860792*x[2] + 0.018736385625077446*x[3] - 0.008392059099094502*x[4]);
    xp[3] = (190.59457176680613*u - 161.57645282794124*x[0] - 0.23872534677018914*x[1] + 1.0884789050298469*x[2] + 0.32394494701618637*x[3] + 0.32518305451736074*x[4]);
    xp[4] = (18.392870361917002*u - 0.43306059549357445*x[0] + 0.60377162139631557*x[1] + 0.62662564832184231*x[2] - 0.48738482327867771*x[3] + 0.98218650191968704*x[4]);

    // Accumulate the output y = C*x + D*u
    y[0] = (182.81664929547824*u - 63.477219815374006*x[0] + 3.5715419988427302*x[1] + 4.1831558072019464*x[2] - 1.0447833362501759*x[3] + 0.27420732436215378*x[4]);

    // Make the predicted state the current state
    for (i=0; i < 5; ++i) {
        x[i] = xp[i];
    }
}
```

## Summary
This tutorial has shown how to follow a workflow that consists of
1. Estimate a process model using experimental data.
2. Design a controller based on the estimated model. We designed a PID controller and one pole-placement controller which was able to cancel the resonances in the system which the PID controllers could not do.
3. Simulate the closed-loop system and analyze its robustness properties. Model uncertainty was considered using the coherence function.
4. Generate C-code for one of the controllers.

Each of these steps is covered in additional detail in the videos available in the playlist [Control systems in Julia](https://youtube.com/playlist?list=PLC0QOsNQS8hZtOQPHdtul3kpQwMOBL8Qc&si=yUrXz5cH4QqTPlR_). See also the tutorial [Control design for a quadruple-tank system](https://help.juliahub.com/juliasimcontrol/dev/examples/quadtank/).