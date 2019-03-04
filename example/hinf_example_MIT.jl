using ControlSystems
using Plots

"""
This is a simple SISO example which was used for debugging the implementation,
as this exact example was use in the lecture notes of the "Principles of Optimal
Control" cours of the MIT OpenCourseWare [1], where the choice of weighting
functions and dynamics gave rise to an H-infinity optimal cotnroller with a
γ of approximately 1.36, where, in our case, we get a controller at 0.93

[1] https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-323-principles-of-optimal-control-spring-2008/lecture-notes/lec15.pdf

The example can be set to visualize and save plots using the two variables
  ShowPlots - true/false (true if plots are to be generated, false for testing)
  filename  - Set to string if files are to be saved, otherwise set a empty list
"""
MakePlots, SavePlots = false, false

# Define the process
G   = tf([200], [0.025,1.0025,10.1,1])

# Sensitivity weight function
M, wB, A = 1.5, 10, 1e-4
WS = tf([1/M, wB],[1, wB*A])

# Output sensitivity weight function
WU = ss(0.1)

# Complementary sensitivity weight function
WT = []

# Form augmented P dynamics in state-space
P = hinfpartition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hinfassumptions(P)

# Synthesize the H-infinity optimal controller
flag, C, gamma = hinfsynthesize(P)

# Extract the transfer functions defining some signals of interest
Pcl, S, CS, T = hinfsignals(P, G, C)

## Plot the specifications
if MakePlots
  specificationplot([S, CS, T], [ss(WS), WU, WT], gamma)
  if SavePlots
    savefig("example_MIT_specifications.pdf")
  end
end

## Plot the closed loop gain from w to z
if MakePlots
  specificationplot(Pcl, gamma; s_labels=["\$\\sigma(P_{cl}(j\\omega))\$"], w_labels=["\$\\gamma\$"])
  if SavePlots
    savefig("example_MIT_clgain.pdf")
  end
end
