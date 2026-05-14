module ControlSystems

using Reexport: Reexport, @reexport
@reexport using ControlSystemsBase
# Explicit imports — only the symbols actually used in this module.
# Derived from `ExplicitImports.print_explicit_imports(ControlSystems)`.
using ControlSystemsBase: ControlSystemsBase, AbstractStateSpace, Continuous,
                          c2d, iscontinuous, ninputs, noutputs, nstates, numeric_type
import DiffEqBase

using LinearAlgebra: LinearAlgebra, diagind, mul!
import OrdinaryDiffEq
import DiffEqCallbacks: SavingCallback, SavedValues
import DelayDiffEq
using SparseArrays: SparseArrays
using StaticArrays: StaticArrays
using RecipesBase: RecipesBase
using Printf: Printf

export Simulator

include("timeresp.jl")
include("simulators.jl")

# The path has to be evaluated upon initial import
const __CONTROLSYSTEMS_SOURCE_DIR__ = dirname(Base.source_path())

end
