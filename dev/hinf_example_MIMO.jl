using ControlSystems
using Plots

"""
This is very dense and horrible example used for debugging. CUrrently doesn't
work due to the absence of enforcing a loopshift in h´the solver - I need to
figure out how to doe this prolperly.
"""

# Define the process dynamics
nx = 30
nu = 8
ny = 8
A = 10*rand(nx,nx)
A -= (abs(maximum(real(eigvals(A))))+0.1) * Matrix{Float64}(I, size(A,1), size(A,2))
B = rand(nx,nu)
C = rand(ny,nx)
D = rand(ny,nu)

G = hInf_bilinear_z2s(hInf_bilinear_s2z(ss(A,B,C,D), 0.005))

# Sensitivity weight function
WSe = tf([0.1,1],[1000,1])
#WS = [WSe 0 0 ; 0 WSe 0 ; 0 0 WSe]
WS = [WSe 0 0 0 ; 0 WSe 0 0 ; 0 0 WSe 0 ; 0 0 0 WSe]
WS = [WS 0*WS; WS*0 WS]

# Output sensitivity weight function
WUe = 0.0001*tf([1,1],[0.1,1])
WU = [WUe 0 0 0 ; 0 WUe 0 0 ; 0 0 WUe 0 ; 0 0 0 WUe]
WU = [WU 0*WU; WU*0 WU]

#WU = [WUe 0 0 ; 0 WUe 0 ; 0 0 WUe ]
# Complementary sensitivity weight function
WTe = tf([10,1],[1,1])
WT = [WTe 0 0 0 ; 0 WTe 0 0 ; 0 0 WTe 0 ; 0 0 0 WTe]
WT = [WT 0*WT; WT*0 WT]

# Form augmented P dynamics in state-space
P = hInf_partition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hInf_assumptions(P)

# Synthesize the H-infinity optimal controller
flag, K, gamma = hInf_synthesize(P; maxIter=20)

Pcl, S, KS, T = hInf_signals(P, G, K)

MakePlots = true
SavePlots = false


fmin = -8
fmax =  8
# Plot the singular values  of the system
f = [10^i for i in range(fmin, stop=fmax, length=10001)]

valPcl  = sigma(Pcl, f)[1];
valS    = sigma(S, f)[1];
valKS   = sigma(KS, f)[1];
valT    = sigma(T, f)[1];

# Visualize the close loop gain
pGain = plot(f, valPcl[:,1], xscale = :log10, yscale = :log10, color = :black, w=2, label="\$\\sigma(P_{w\\rightarrow z}(j\\omega))\$")
plot!(f, gamma*ones(size(f)), xscale = :log10, yscale = :log10, color = :black, w=3, style=:dot, label="\$\\gamma\$")
if size(valPcl,2) > 1
  for ii = 2:size(valPcl,2)
    plot!(f, valPcl[:,ii], xscale = :log10, yscale = :log10, color = :black, w=1)
  end
end

# Plot the sensitivity functions of interest
pSigma = plot(f, valS[:,1], xscale = :log10, yscale = :log10, color = :red, w=3, label="\$\\sigma(S(j\\omega))\$")
if size(valS,2) > 1
  for ii = 2:size(valS,2)
    plot!(f, valS[:,ii], xscale = :log10, yscale = :log10, color = :red, w=1, label="")
  end
end
plot!(f, valKS[:,1], xscale = :log10, yscale = :log10, color = :green, w=3, label="\$\\sigma(C(j\\omega)S(j\\omega)))\$")
if size(valKS,2) > 1
  for ii = 2:size(valKS,2)
    plot!(f, valKS[:,ii], xscale = :log10, yscale = :log10, color = :green, w=1, label="")
  end
end
plot!(f, valT[:,1], xscale = :log10, yscale = :log10, color = :blue, w=3, label="\$\\sigma(T(j\\omega)))\$")
if size(valT,2) > 1
  for ii = 2:size(valT,2)
    plot!(f, valT[:,ii], xscale = :log10, yscale = :log10, color = :blue, w=1, label="")
  end
end

# Visualize the weighting functions
if isa(WS, LTISystem) || isa(WS, Number)
  valiWS = sigma(gamma/WS[1,1], f)[1]
  plot!(f, valiWS, xscale = :log10, yscale = :log10, color = :red, w=2, style=:dot, label="\$\\gamma\\sigma(W_S(j\\omega)^{-1})\$")
end
if isa(WU, LTISystem) || isa(WU, Number)
  valiWU = sigma(gamma/WU[1,1], f)[1]
  plot!(f, valiWU, xscale = :log10, yscale = :log10, color = :green, w=2, style=:dot, label="\$\\gamma\\sigma(W_U(j\\omega)^{-1})\$")
end
if isa(WT, LTISystem) || isa(WT, Number)
  valiWT = sigma(gamma/WT[1,1], f)[1]
  plot!(f, valiWT, xscale = :log10, yscale = :log10, color = :blue, w=2, style=:dot, label="\$\\gamma\\sigma(W_T(j\\omega)^{-1})\$")
end

l = @layout [ a ; b ]
plt=plot(pGain, pSigma, layout=l, size=(600,600))
gui()
