# % A flexible transmission
# % A not so succesful design
# format compact

# % Input your RST controller here
using ControlSystemsBase, DSP, Plots, Statistics, MonteCarloMeasurements, RobustAndOptimalControl
unsafe_comparisons(true)
r = [1, -1]
s = 0.02*[1, 0]
t = [1, 0 ]
t *= sum(s)/sum(t)


fs = 20
h = 1/fs

fv  =  0.01:0.01:10
wv = 2 .* pi .* fv
b0 = [0, 0, 0, 0.28261, 0.50666]
a0 = [1, -1.41833, 1.58939, -1.31608, 0.88642]
b05 = [0, 0, 0, 0.10276, 0.18123]
a05 = [1, -1.99185, 2.20265, -1.84083, 0.89413]
b1 = [0, 0, 0, 0.06408, 0.10407]
a1 = [1, -2.09679, 2.31962, -1.93353, 0.87129]
P0 = tf(b0,a0,h) |> ss
P05 = tf(b05,a05,h) |> ss
P1 = tf(b1,a1,h) |> ss

P = ss2particles([P0, P05, P1])



function plotTF_uncertain(G)
    # Plot frequency response for uncertain system G
    C = abs.(freqrespv(G, wv))
    plot(fv, 20log10.(C))
    return C
end

plotTF_uncertain(P)
title!("The Process")
xlabel!("Frequency Hz")
ylabel!("dB")


Cfb = tf(s,r,h)
Cff = tf(t,r,h)

# % Specification A and B
T = 0:h:5
y_particles = step(Cff*feedback(P, Cfb), T).y[:]
y_matrix = Matrix(y_particles)  # Convert to matrix: particles × time
plot(T, y_matrix')
plot!([0, 5],[1.1, 1.1],c=:black,linewidth=2)
plot!([1.15, 5],[0.9, 0.9],c=:black,linewidth=2)
plot!([1.15, 5],[1.1, 1.1],c=:black,linewidth=2)
title!("Step Response")

# Calculate worst-case metrics across all particles
risetime = [maximum(T[abs.(y_matrix[i,:].-1) .> 0.1])-0.15 for i in 1:size(y_matrix,1)]  # Compensation for 3*h
risetimemax = maximum(risetime)
PASSA = risetimemax <= 1

peaks = [maximum(y_matrix[i,:]) for i in 1:size(y_matrix,1)]
peakmax = maximum(peaks)
PASSB = (peakmax <= 1.1)


# % Specification C - Disturbance response
# Use feedback to get the disturbance transfer function directly
sysd = feedback(P, Cfb)
d_particles = step(sysd, T).y[:]
d_matrix = Matrix(d_particles)  # Convert to matrix: particles × time

# Calculate worst-case metrics
peakd = [maximum(abs, d_matrix[i, :]) for i in axes(d_matrix, 1)]
timed = [maximum(T[abs.(d_matrix[i, :]) .> 0.1*peakd[i]]) for i in axes(d_matrix, 1)]
timedmax = maximum(timed)
PASSC = timedmax < 1.2

plot(T, d_matrix')
peakd_max = maximum(peakd)
plot!([1.2, maximum(T)], peakd_max*[0.1, 0.1], c=:black, linewidth=2, legend=false)
plot!([1.2, maximum(T)], -peakd_max*[0.1, 0.1], c=:black, linewidth=2, legend=false)
title!("1/A output disturbance Response")


# % Specification D
PASSD = abs(evalfr(Cfb,1)[1]) > 1e5

# % Specification E and F - Sensitivity function
Sens = feedback(1, P*Cfb)

S_particles = abs.(freqrespv(Sens, wv))
S_matrix = Matrix(S_particles)  # Convert to matrix: particles × frequency
plot(fv, 20log10.(S_matrix'))
plot!([0, 0.2, 0.2, 10],[0, 0, 6, 6],c=:black,linewidth=2)
title!("Sensitivity function")

# Calculate worst-case metrics across all particles
freqe = [minimum(fv[S_matrix[i, :] .> 1]) for i in axes(S_matrix, 1)]
freqemin = minimum(freqe)
PASSE = freqemin > 0.2

peake = [20log10(maximum(S_matrix[i, :])) for i in axes(S_matrix, 1)]
peakemax = maximum(peake)
PASSF = peakemax < 6

# % Specification G - Stability margins
# For uncertain systems, calculate margins for each particle
Margins = []
DM = []

for i in axes(P.A, 3)
    Pi = sys_from_particles(P, i)
    push!(Margins, margin(Pi*Cfb, allMargins=false))
    push!(DM, delaymargin(Pi*Cfb))
end

DMmin = minimum(DM)
PASSG = DMmin > 0.8

# % Specification H - Control effort sensitivity
CSens = feedback(Cfb, P)

CS_particles = abs.(freqrespv(CSens, wv))
CS_matrix = Matrix(CS_particles)  # Convert to matrix: particles × frequency
plot(fv, 20log10.(CS_matrix'), xscale=:log10)
plot!([8, 10],[10, 10],c=:black,linewidth=2)
title!("CS")

ind8 = (wv .> 2pi*8)
peakh = [20log10(maximum(CS_matrix[i, ind8])) for i in axes(CS_matrix, 1)]
peakhmax = maximum(peakh)
PASSH = peakhmax < 10


PASSvector =[PASSA PASSB PASSC PASSD PASSE PASSF PASSG PASSH][:]
# % Check stability for all particles
stable(m) = m[2][1] > 0 && m[4][1] > 0
Allstable = all(stable(m) for m in Margins)


Scorev = round(100*min(1,max(0,max(2-risetimemax, (1.2-peakmax)/0.1, (1.4-timedmax)/0.2, PASSD, (freqemin-0.1)/0.1, (12-peakemax)/6, (DMmin-0.4)/0.4, (20-peakhmax)/10))))
SCORE = Allstable*Scorev