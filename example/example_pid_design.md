
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
