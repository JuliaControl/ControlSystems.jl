using ControlSystems
using Plots
using LinearAlgebra
"""
This is a simple SISO example with integrator dynamics corresponding to the
quad tank process in the lab.

The example can be set to visualize and save plots using the variables
  MakePlots - true/false (true if plots are to be generated, false for testing)
  SavePlots - true/false (true if plots are to be saved, false for testing)
"""
MakePlots = true
SavePlots = true

# Define the proces parameters
k1, k2, kc, g = 3.33, 3.35, 0.5, 981
A1, A3, A2, A4 = 28, 28, 32, 32
a1, a3, a2, a4= 0.071, 0.071, 0.057, 0.057
h01, h02, h03, h04 = 12.4, 12.7, 1.8, 1.4
T1, T2 = (A1/a1)*sqrt(2*h01/g), (A2/a2)*sqrt(2*h02/g)
T3, T4 = (A3/a3)*sqrt(2*h03/g), (A4/a4)*sqrt(2*h04/g)
c1, c2 = (T1*k1*kc/A1), (T2*k2*kc/A2)
gamma1, gamma2 = 0.7, 0.6

# Define the process dynamics
A = [-1/T1     0 A3/(A1*T3)          0;
     0     -1/T2          0 A4/(A2*T4);
     0         0      -1/T3          0;
     0         0          0      -1/T4];
B = [gamma1*k1/A1     0;
     0                gamma2*k2/A2;
     0                (1-gamma2)*k2/A3;
     (1-gamma1)*k1/A4 0              ];
C = [kc 0 0 0;
     0 kc 0 0];
D = zeros(2,2)
G = ss(A,B,C,D);

# Sensitivity weight function
WSelement = 100*tf([1,1],[1000,1])
WS = [WSelement 0; 0 WSelement]
iWSelement = 1/WSelement
iWS = [iWSelement 0; 0 iWSelement]

# Output sensitivity weight function
WUelement = 5*tf([1,1],[0.1,1]) ##
WUelement = ss(0.01)
WU = [WUelement 0; 0 WUelement]
iWUelement = 1/WUelement
iWU = [iWUelement 0; 0 iWUelement]

# Complementary sensitivity weight function
WTelement = tf([10,0.1],[1,1])
WT  = [WTelement 0; 0 WTelement]
iWTelement = 1/WTelement
iWT = [iWTelement 0; 0 iWTelement]

# Form augmented P dynamics in state-space
P = hInf_partition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hInf_assumptions(P)

# Synthesize the H-infinity optimal controller
flag, C, gamma = hInf_synthesize(P)

Pcl, S, CS, T = hInf_signals(P, G, C)

if MakePlots
  # Visualize
  frequencies = [10^i for i in range(-6, stop=6, length=1001)];

  pGain = plot(frequencies, ones(size(frequencies))*gamma, scale = :log10, yscale = :log10, color = :black, w=2, style=:dot, label="\$\\gamma\$")
  valPcl  = sigma(Pcl, frequencies)[1];
  plot!(frequencies, valPcl[:,1], scale = :log10, yscale = :log10, color = :black, w=1, label="\$\\sigma(P_{w\\rightarrow z}(j\\omega))\$")
  plot!(frequencies, valPcl[:,2], scale = :log10, yscale = :log10, color = :black, w=1, label="")

  valT  = sigma(T, frequencies)[1];
  pSigma = plot(frequencies, valT[:,1], scale = :log10, yscale = :log10, color = :black, w=1, label="\$\\sigma(T(j\\omega))\$")
  plot!(frequencies, valT[:,2], scale = :log10, yscale = :log10, color = :black, w=1, label="")
  valCS  = sigma(CS, frequencies)[1];
  plot!(frequencies, valCS[:,1], scale = :log10, yscale = :log10, color = :blue, w=1, label="\$\\sigma(C(j\\omega)S(j\\omega))\$")
  plot!(frequencies, valCS[:,2], scale = :log10, yscale = :log10, color = :blue, w=1, label="")
  valS = sigma(S, frequencies)[1];
  plot!(frequencies, valS[:,1], scale = :log10, yscale = :log10, color = :red, w=1, label="\$\\sigma(S(j\\omega))\$")
  plot!(frequencies, valS[:,2], scale = :log10, yscale = :log10, color = :red, w=1, label="")

  valPiWT  = sigma(gamma*iWT, frequencies)[1];
  plot!(frequencies, valPiWT[:,1], scale = :log10, yscale = :log10, color = :black, w=2, style=:dot, label="\$\\gamma\\sigma(W_T^{-1}(j\\omega))\$")#
  valPiWU  = sigma(gamma*iWU, frequencies)[1];
  plot!(frequencies, valPiWU[:,1], scale = :log10, yscale = :log10, color = :blue, w=2, style=:dot, label="\$\\gamma\\sigma(W_U^{-1}(j\\omega))\$")
  valPiWS  = sigma(gamma*iWS, frequencies)[1];
  plot!(frequencies, valPiWS[:,1], scale = :log10, yscale = :log10, color = :red, w=2, style=:dot, label="\$\\gamma\\sigma(W_S^{-1}(j\\omega))\$")

  times = [i for i in range(0, stop=300, length=10000)]

  stepy, stept, stepx = step(T, times)

  pStep1=plot(stept, stepy[:,1,1], color = :blue, w=2, label="\$u_1\\rightarrow y_1\$")
  pStep2=plot(stept, stepy[:,1,2], color = :blue, w=2, label="\$u_1\\rightarrow y_2\$", ylims = (-0.5,1.1))
  pStep3=plot(stept, stepy[:,2,1], color = :blue, w=2, label="\$u_2\\rightarrow y_1\$", ylims = (-0.5,1.1))
  pStep4=plot(stept, stepy[:,2,2], color = :blue, w=2, label="\$u_2\\rightarrow y_2\$")

  if SavePlots
    plot(pGain, size=(600, 300), legend=:right)
    savefig("gainPlotQuadTank.pdf")

    plot(pSigma, size=(600, 300), legend=:left)
    savefig("sigmaPlotQuadTank.pdf")

    l = @layout [ c d e f ]
    plt=plot(pStep1, pStep2, pStep3, pStep4, layout=l, size=(1000,250))
    savefig("stepPlotQuadTank.pdf")
  else
    l = @layout [ a ; b ; c d ; e f ]
    plt=plot(pGain, pSigma, pStep1, pStep2, pStep3, pStep4, layout=l, size=(600,600))
  end
end
