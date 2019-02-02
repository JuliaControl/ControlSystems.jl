using ControlSystems, Plots
gr()

sys = feedback(1.0, ss(-1.0, 2, 1, 0) * (delay(2.0) + delay(3.0) + delay(2.5)))
sys = feedback(ss(-1.0, 1, 1, 0), delay(1.0))

t = 0:0.02:8
@time y, t, x = lsim(sys, t-> [t>=0 ? 1.0 : 0.0], t)
@time y, t, x = lsim(sys, [1.0], t)
@time y, t, x = lsim(sys, (out, t) -> (out .= (t>=0 ? 1.0 : 0.0)), t)
@time y, t, x = lsim(sys, (out, t) -> (out[1] = (t>=0 ? 1.0 : 0.0)), t)

function u0(out,t)
    if t > 0
        out[1] = 1
    else
        out[1] = 0
    end
    return
end

@time y, t, x = lsim(sys, u0, t)

plot(t, y')

s = tf("s")
P = delay(2.6)*ss((s+3.0)/(s^2+0.3*s+1))
C = 0.06 * ss(1.0 + 1/s);
P*C
T = feedback(P*C,1.0)

t = 0:0.1:70
y, t, x = lsim(T, t -> (t<0 ? 0 : 1 ), t)
plot(t, y, c = :blue)

w = 10 .^ (-2:0.01:2)
marginplot(P*C, w)
marginplot(P*C)

notch = ss(tf([1, 0.2, 1],[1, .8, 1]));
C = ss(0.05 * (1 + 1/s));
Tnotch = feedback(P*C*notch, 1.0)

stepplot(Tnotch)

y, t, x = step(C, method=:zoh)

y2, t2, x2 = step(Tnotch)
stepplot(Tnotch)

stepplot(Tnotch, 40, 0.1)

stepplot(T, 100)

G = delay(5)/(s+1)
T = feedback(G, 0.5)
w = 10 .^ (-2:0.01:3)
bodeplot(T, w, plotphase=false)

# Test conversion, promotion
delay(1,Int64) + 3.5

G = 1 + 0.5 * delay(3)
w = 10 .^(-2:0.001:2)
bodeplot(G, w, plotphase=false)

G = delay(1) * ((0.8*s^2+s+2)/(s^2+s))
T = feedback(G,1)
# Not possible with direct term
stepplot(T)

bodeplot(T)

G = 1/(s+1) + delay(4)
T = feedback(1,G)
# Not possible to lsim with direct term
stepplot(T)
bodeplot(T)

s = tf("s")
