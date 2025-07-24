module ControlSystems

using Reexport
@reexport using ControlSystemsBase
using ControlSystemsBase: issiso, ninputs, noutputs, nstates, numeric_type


using LinearAlgebra
import OrdinaryDiffEq
import LinearAlgebra: BlasFloat
import Hungarian # For pole assignment in rlocusplot
import DiffEqCallbacks: SavingCallback, SavedValues
import DelayDiffEq
using SparseArrays
using StaticArrays
using RecipesBase
using Printf

export Simulator

include("timeresp.jl")
include("simulators.jl")
include("root_locus.jl")

# The path has to be evaluated upon initial import
const __CONTROLSYSTEMS_SOURCE_DIR__ = dirname(Base.source_path())

end
