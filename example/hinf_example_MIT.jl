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
MakePlots = false
SavePlots = false

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
P = hInf_partition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hInf_assumptions(P)

# Synthesize the H-infinity optimal controller
flag, C, gamma = hInf_synthesize(P)

# Extract the transfer functions defining some signals of interest
Pcl, S, CS, T = hInf_signals(P, G, C)

# TODO remove hack for visualizing plots, should be made into some kind of recepie
if MakePlots
  include("hinf_utilities.jl")
  if SavePlots; filename = "example_MIT.pdf"; else; filename=[]; end
  visualize_synthesis(Pcl, S, CS, T, gamma, tmax=1)
end
