function makePlots()

println("Generating plots")

plotsDir = (pwd()[end-3:end] == "docs") ? "build/plots" : "docs/build/plots"
mkdir(plotsDir)
Plots.pyplot()


# PID design functions
P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
Plots.savefig(plotsDir*"/pidgofplot.svg")

ωp = 0.8
kp,ki,C = loopshapingPI(P,ωp,phasemargin=60, doplot=false)
gangoffourplot(P, [tf(1), C])
Plots.savefig(plotsDir*"/pidgofplot2.svg")
nyquistplot([P, P*C])
Plots.savefig(plotsDir*"/pidnyquistplot.svg")


ωp = 2
kp,ki,C60 = loopshapingPI(P,ωp,rl=1,phasemargin=60, doplot=true)
gangoffourplot(P, [tf(1), C60])
Plots.savefig(plotsDir*"/pidgofplot3.svg")
nyquistplot([P, P*C60])
Plots.savefig(plotsDir*"/pidnyquistplot2.svg")

# Advanced pole placement
ζ = 0.2
ω = 1
B = [1]
A   = [1, 2ζ*ω, ω^2]
P  = tf(B,A)
# Control design
ζ0 = 0.7
ω0 = 2
Am = [1, 2ζ0*ω0, ω0^2]
Ao = conv(2Am, [1/2, 1]) # Observer polynomial
AR = [1,0] # Force the controller to contain an integrator

B⁺  = [1] # The process numerator polynomial can be facored as B = B⁺B⁻ where B⁻ contains the zeros we do not want to cancel (non-minimum phase and poorly damped zeros)
B⁻  = [1]
Bm  = conv(B⁺, B⁻) # In this case, keep the entire numerator polynomial of the process

R,S,T = rstc(B⁺,B⁻,A,Bm,Am,Ao,AR) # Calculate the 2-DOF controller polynomials

Gcl = tf(conv(B,T),zpconv(A,R,B,S)) # Form the closed loop polynomial from reference to output

stepplot([P,Gcl]) # Visualize the open and closed loop responses.
Plots.savefig(plotsDir*"/ppstepplot.svg")
gangoffourplot(P, tf(-S,R)) # Plot the gang of four to check that all tranfer functions are OK
Plots.savefig(plotsDir*"/ppgofplot.svg")

P1 = "exp(-sqrt(s))"
f1 = stabregionPID(P1,logspace(-5,1,1000)); Plots.savefig(plotsDir*"/stab1.svg")
P2 = "100*(s+6).^2./(s.*(s+1).^2.*(s+50).^2)"
f2 = stabregionPID(P2,logspace(-5,2,1000)); Plots.savefig(plotsDir*"/stab2.svg")
P3 = tf(1,[1,1])^4
f3 = stabregionPID(P3,logspace(-5,0,1000)); Plots.savefig(plotsDir*"/stab3.svg")




# PID plots

P = tf([1.],[1., 1])
ζ = 0.5 # Desired damping
ws = logspace(-1,2,8) # A vector of closed-loop bandwidths
kp = 2*ζ*ws-1 # Simple pole placement with PI given the closed-loop bandwidth, the poles are placed in a butterworth pattern
ki = ws.^2
pidplots(P,:nyquist,;kps=kp,kis=ki, ω= logspace(-2,2,500))
Plots.savefig(plotsDir*"/pidplotsnyquist1.svg")
pidplots(P,:gof,;kps=kp,kis=ki, ω= logspace(-2,2,500))
Plots.savefig(plotsDir*"/pidplotgof1.svg")

kp = linspace(-1,1,8) # Now try a different strategy, where we have specified a gain crossover frequency of 0.1 rad/s
ki = sqrt(1-kp.^2)/10
pidplots(P,:nyquist,;kps=kp,kis=ki)
Plots.savefig(plotsDir*"/pidplotsnyquist2.svg")
pidplots(P,:gof,;kps=kp,kis=ki)
Plots.savefig(plotsDir*"/pidplotsgof2.svg")

#AR = [1,0]
#B⁺  = [1]
#B⁻  = [1]
#Bm  = conv(B⁺, B⁻)
#R,S,T = rstc(B⁺,B⁻,A,bm,Am,Ao,AR)
#Gcl = tf(conv(B,T),zpconv(A,R,B,S))
#f1 = stepplot([P,Gcl])
#f2 = gangoffourplot(P, tf(-S,R))
#Plots.savefig(f1, "$plotsDir/rstcstepplot.svg")
#Plots.savefig(f2, "$plotsDir/rstcgofplot.svg")
end
