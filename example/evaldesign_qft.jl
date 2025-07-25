# % A flexible transmission
# % A not so succesful design
# format compact

# % Input your RST controller here
using ControlSystemsBase, DSP, Plots
r = [1, -1]
s = 0.02*[1, 0]
t = [1, 0 ]
t *= sum(s)/sum(t)


fs = 20
h = 1/fs

fv  =  0.01:0.01:10
wv = 2Ï€*fv
b0 = [0, 0, 0, 0.28261, 0.50666]
a0 = [1, -1.41833, 1.58939, -1.31608, 0.88642]
b05 = [0, 0, 0, 0.10276, 0.18123]
a05 = [1, -1.99185, 2.20265, -1.84083, 0.89413]
b1 = [0, 0, 0, 0.06408, 0.10407]
a1 = [1, -2.09679, 2.31962, -1.93353, 0.87129]
P0 = tf(b0,a0,h)
P05 = tf(b05,a05,h)
P1 = tf(b1,a1,h)

function plotTFs(G0, G05, G1)
    C0  = abs(squeeze(freqresp(G0,wv)[1],(1,2)))
    C05 = abs(squeeze(freqresp(G05,wv)[1],(1,2)))
    C1  = abs(squeeze(freqresp(G1,wv)[1],(1,2)))
    plot(fv,20log10(C0))
    plot!(fv,20log10(C05),c=:red)
    plot!(fv,20log10(C1),c=:black)
    return C0, C05, C1
end

plotTFs(tf(b0,a0,1/fs),tf(b05,a05,1/fs),tf(b1,a1,1/fs))
title!("The Process")
xlabel!("Frequency Hz")
ylabel!("dB")


Cfb = tf(s,r,h)
Cff = tf(t,r,h)

# % Specification A and B
T = 0:h:5
y0 = step(Cff*feedback(P0, Cfb),T)[1]
y05 = step(Cff*feedback(P05, Cfb),T)[1]
y1 = step(Cff*feedback(P1, Cfb),T)[1]
plot(T,y0)
plot!(T,y05,c=:red)
plot!(T,y1,c=:black)
plot!([0, 5],[1.1, 1.1],c=:black,linewidth=2)
plot!([1.15, 5],[0.9, 0.9],c=:black,linewidth=2)
plot!([1.15, 5],[1.1, 1.1],c=:black,linewidth=2)
title!("Step Response")

risetime0 = maximum(T[find(abs(y0-1).>0.1)])-0.15     # Compensation for 3*h
risetime05 = maximum(T[find(abs(y05-1).>0.1)])-0.15
risetime1= maximum(T[find(abs(y1-1).>0.1)])-0.15
risetimemax =  maximum([risetime0,risetime05,risetime1])
PASSA =  maximum([risetime0,risetime05,risetime1])<=1

peak0 = maximum(y0)
peak05 = maximum(y05)
peak1 = maximum(y1)
peakmax = maximum([peak0, peak05, peak1])
PASSB = (peakmax <= 1.1)


# % Specification C
sysd0 = tf(conv(r,[1, 0, 0, 0, 0]),conv(a0,r)+conv(b0,s),h)
sysd05 = tf(conv(r,[1, 0, 0, 0, 0]),conv(a05,r)+conv(b05,s),h)
sysd1 = tf(conv(r,[1, 0, 0, 0, 0]),conv(a1,r)+conv(b1,s),h)
d0 = step(sysd0,T)[1]
d05 = step(sysd05,T)[1]
d1 = step(sysd1,T)[1]
peakd0 = maximum(abs(d0))
peakd05 = maximum(abs(d05))
peakd1 = maximum(abs(d1))
timed0 = maximum(T[find(abs(d0).>0.1*peakd0)])
timed05 = maximum(T[find(abs(d05).>0.1*peakd05)])
timed1= maximum(T[find(abs(d1).>0.1*peakd1)])
timedmax = maximum([timed0,timed05,timed1])
PASSC = timedmax < 1.2
plot(T,[d0 d05 d1],c=[:blue,:red,:black]',lab=[0,0.5,1]')
plot!([1.2, maximum(T)],peakd0*[0.1, 0.1],c=:blue,linewidth=2,legend=false)
plot!([1.2, maximum(T)],-peakd0*[0.1, 0.1],c=:blue,linewidth=2,legend=false)
plot!([1.2, maximum(T)],peakd05*[0.1, 0.1],c=:red,linewidth=2,legend=false)
plot!([1.2, maximum(T)],-peakd05*[0.1, 0.1],c=:red,linewidth=2,legend=false)
plot!([1.2, maximum(T)],peakd1*[0.1, 0.1],c=:black,linewidth=2,legend=false)
plot!([1.2, maximum(T)],-peakd1*[0.1, 0.1],c=:black,linewidth=2,legend=false)
title!("1/A output disturance Response")


# % Specification D
PASSD = abs(evalfr(Cfb,1))[1] > 1e5

# % Specification E and F

Sens0  = feedback(1, P0*Cfb)
Sens05 = feedback(1, P05*Cfb)
Sens1  = feedback(1, P1*Cfb)


S0, S05, S1 =plotTFs(Sens0, Sens05, Sens1)
# axis([0, 10, -10, 10])
plot!([0, 0.2, 0.2, 10],[0, 0, 6, 6],c=:black,linewidth=2)
title!("Sensitivity function")
plot!(xscale=:log10)

freqe0 = minimum(fv[find(S0.>1)])
freqe05 = minimum(fv[find(S05.>1)])
freqe1 = minimum(fv[find(S1.>1)])
freqemin =  minimum([freqe0 freqe05 freqe1])
PASSE = freqemin > 0.2

peake0 = 20log10(maximum(S0))
peake05 = 20log10(maximum(S05))
peake1 = 20log10(maximum(S1))

peakemax = maximum([peake0 peake05 peake1])
PASSF = peakemax < 6

# % Specification G

function delaymargin(G)
    # Phase margin in radians divided by cross-over frequency in rad/s.
    m   = margin(G,allMargins=true)
    Ï•â‚˜,i= findmin(m[4])
    Ï•â‚˜ *= Ï€/180
    Ï‰Ï•â‚˜ = m[3][i]
    dâ‚˜  = (Ï•â‚˜/Ï‰Ï•â‚˜/h)[1] # Give delay margin in number of sample times, as matlab does
    dâ‚˜ <= 0 && warn("dâ‚˜ <= 0")
    dâ‚˜
end

Margins0 = margin(P0*Cfb,allMargins=false)
Margins05 = margin(P05*Cfb,allMargins=false)
Margins1 = margin(P1*Cfb,allMargins=false)

DM0 = delaymargin(P0*Cfb)
DM05 = delaymargin(P05*Cfb)
DM1 = delaymargin(P1*Cfb)

DMmin = minimum([DM0 DM05 DM1])
PASSG = DMmin>0.8

# % Specification H

CSens0  = Cfb/(1+P0*Cfb)
CSens05 = Cfb/(1+P05*Cfb)
CSens1  = Cfb/(1+P1*Cfb)
CS0, CS05, CS1 = plotTFs(CSens0,CSens05,CSens1)
plot!([8, 10],[10, 10],c=:black,linewidth=2, xscale=:log10)
# axis([0, 10, -15, 15])
title!("CS")

ind8 = find(wv.>2pi*8)
peakh0 = 20log10(maximum(CS0[ind8]))
peakh05 = 20log10(maximum(CS05[ind8]))
peakh1 = 20log10(maximum(CS1[ind8]))
peakhmax = maximum([peakh0, peakh05, peakh1])
PASSH = peakhmax < 10


PASSvector =[PASSA PASSB PASSC PASSD PASSE PASSF PASSG PASSH][:]
# % And lets check stability also !
stable(m) = m[2][1] > 0 && m[4][1] > 0
Allstable = stable(Margins0) && stable(Margins05) && stable(Margins1)


Scorev = round(100*min(1,max(0,[2-risetimemax  (1.2-peakmax)/0.1  (1.4-timedmax)/0.2  PASSD (freqemin-0.1)/0.1 (12-peakemax)/6 (DMmin-0.4)/0.4 (20-peakhmax)/10]  )))
SCORE = Allstable*mean(Scorev)