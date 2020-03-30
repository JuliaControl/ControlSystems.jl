# PID design functions
By plotting the gang of four under unit feedback for the process
```julia
P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
```
we notice that the sensitivity function is a bit too high at around frequencies ω = 0.8 rad/s. Since we want to control the process using a simple PI-controller, we utilize the
function `loopshapingPI` and tell it that we want 60 degrees phase margin at this frequency. The resulting gang of four is plotted for both the constructed controller and for unit feedback.

```julia
ωp = 0.8
kp,ki,C = loopshapingPI(P,ωp,phasemargin=60, doplot=true)
```

We could also cosider a situation where we want to create a closed-loop system with the bandwidth ω = 2 rad/s, in which case we would write something like
```julia
ωp = 2
kp,ki,C60 = loopshapingPI(P,ωp,rl=1,phasemargin=60, doplot=true)
```
Here we specify that we want the Nyquist curve `L(iω) = P(iω)C(iω)` to pass the point `|L(iω)| = rl = 1,  arg(L(iω)) = -180 + phasemargin = -180 + 60`
The gang of four tells us that we can indeed get a very robust and fast controller with this design method, but it will cost us significant control action to double the bandwidth of all four poles. 

# Advanced pole-zero placement
This example illustrates how we can perform advanced pole-zero placement. The task is to make the process a bit faster and damp the poorly damped poles.

Define the process
```julia
ζ = 0.2
ω = 1

B = [1]
A   = [1, 2ζ*ω, ω^2]
P  = tf(B,A)
```

Define the desired closed loop response, calculate the controller polynomials and simulate the closed-loop system. The design utilizes an observer poles twice as fast as the closed-loop poles. An additional observer pole is added in order to get a casual controller when an integrator is added to the controller. 
```julia
# Control design
ζ0 = 0.7
ω0 = 2
Am = [1, 2ζ0*ω0, ω0^2]
Ao = conv(2Am, [1/2, 1]) # Observer polynomial, add extra pole due to the integrator
AR = [1,0] # Force the controller to contain an integrator ( 1/(s+0) )

B⁺  = [1] # The process numerator polynomial can be facored as B = B⁺B⁻ where B⁻ contains the zeros we do not want to cancel (non-minimum phase and poorly damped zeros)
B⁻  = [1]
Bm  = conv(B⁺, B⁻) # In this case, keep the entire numerator polynomial of the process

R,S,T = rstc(B⁺,B⁻,A,Bm,Am,Ao,AR) # Calculate the 2-DOF controller polynomials

Gcl = tf(conv(B,T),zpconv(A,R,B,S)) # Form the closed loop polynomial from reference to output, the closed-loop characteristic polynomial is AR + BS, the function zpconv takes care of the polynomial multiplication and makes sure the coefficient vectores are of equal length

stepplot([P,Gcl]) # Visualize the open and closed loop responses.
gangoffourplot(P, tf(-S,R)) # Plot the gang of four to check that all tranfer functions are OK
```



# Stability boundary for PID controllers
The stability boundary, where the transfer function `P(s)C(s) = -1`, can be plotted with the command `stabregionPID`. The process can be given in string form or as a regular LTIsystem.

```julia
P1 = "exp(-sqrt(s))"
f1 = stabregionPID(P1,exp10.(range(-5, stop=1, length=1000)))
P2 = "100*(s+6).^2./(s.*(s+1).^2.*(s+50).^2)"
f2 = stabregionPID(P2,exp10.(range(-5, stop=2, length=1000)))
P3 = tf(1,[1,1])^4
f3 = stabregionPID(P3,exp10.(range(-5, stop=0, length=1000)))
```




# PID plots
This example utilizes the function `pidplots`, which accepts vectors of PID-parameters and produces relevant plots. The task is to take a system with bandwidth 1 rad/s and produce a closed-loop system with bandwidth 0.1 rad/s. If one is not careful and proceed with pole placement, one easily get a system with very poor robustness.
```julia
P = tf([1.],[1., 1])

ζ = 0.5 # Desired damping

ws = logspace(-1,2,8) # A vector of closed-loop bandwidths
kp = 2*ζ*ws-1 # Simple pole placement with PI given the closed-loop bandwidth, the poles are placed in a butterworth pattern
ki = ws.^2
pidplots(P,:nyquist,:gof;kps=kp,kis=ki, ω= logspace(-2,2,500)) # Request Nyquist and Gang-of-four plots (more plots are available, see ?pidplots )

kp = linspace(-1,1,8) # Now try a different strategy, where we have specified a gain crossover frequency of 0.1 rad/s
ki = sqrt(1-kp.^2)/10
pidplots(P,:nyquist,:gof;kps=kp,kis=ki)
```
