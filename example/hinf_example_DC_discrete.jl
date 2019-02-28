using Plots
using ControlSystems
"""
This is a simple SISO example with a pole in the origin, corresponding to the
DC servos used in the Lund laboratories. It serves to exeplify how the syntheis
can be done for simple SISO systems, and also demonstrates how we chan verify
if the problem is feasible to solve using the ARE method.

The example can be set to visualize and save plots using the variables
  MakePlots - true/false (true if plots are to be generated, false for testing)
  SavePlots - true/false (true if plots are to be saved, false for testing)
"""
MakePlots = false
SavePlots = false

# Define the process
ts = 0.01
epsilon = 1e-5
Gd = ss(c2d(tf([11.2], [1, 0.12]) * tf([1], [1, epsilon]), ts))

# Sensitivity weight function
M, wB, A = 1.5, 20.0, 1e-8
WS = tf([1/M, wB],[1, wB*A])

# Output sensitivity weight function
WU = ss(1)

# Complementary sensitivity weight function
WT = []

# Create continuous time approximation of the process
hInf_bilinear_z2s(ss(Gd))

# Form the P in the LFT Fl(P,C) as a partitioned state-space object
Pc = hInf_partition(Gc, WS, WU, WT)

# Check if the problem is feasible
flag = hInf_assumptions(Pc)

# Synthesize the H-infinity optimal controller
flag, Cc, gamma = hInf_synthesize(Pc)

# Extract the transfer functions defining some signals of interest, but do so
# using discrete equivalent of the continuous time objects Pc, Cc and Gc
PclD, SD, CSD, TD = hInf_signals(
  hInf_bilinear_s2z(Pc, ts),
  hInf_bilinear_s2z(Gc, ts),
  hInf_bilinear_s2z(Cc, ts)
)

# This solution is a bit hacky and should be revised
Pcl = ss(PclD.A, PclD.B, PclD.C, PclD.D, ts)
S   = ss(SD.A, SD.B, SD.C, SD.D, ts)
CS  = ss(CSD.A, CSD.B, CSD.C, CSD.D, ts)
T   = ss(TD.A, TD.B, TD.C, TD.D, ts)

# TODO remove hack for visualizing plots, should be made into some kind of recepie
if MakePlots
  include("hinf_utilities.jl")
  if SavePlots; filename = "example_DC.pdf"; else; filename=[]; end
  visualize_synthesis(Pcl, S, CS, T, gamma; filename=filename, tmax=4)
end
