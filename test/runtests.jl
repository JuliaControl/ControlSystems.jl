using ControlSystems
using Test, LinearAlgebra, Random
using Aqua
# @testset "Aqua" begin
#     Aqua.test_all(ControlSystems;
#         ambiguities = false, # casues 100s of hits in all dependencies
#         stale_deps = true,  # Aqua complains about itself https://github.com/JuliaTesting/Aqua.jl/issues/78 
#     )
# end

# Helper to call dev on one of the packages in folder /lib
function dev_subpkg(subpkg)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", subpkg)
    Pkg.develop(PackageSpec(path = subpkg_path))
end

dev_subpkg("ControlSystemsBase") # Always dev this package to test with the latest code

const GROUP = ENV["GROUP"] # Get the GROUP attribute from the test.yml file

if GROUP == "ControlSystems"
    include("runtests_controlsystems.jl")
else
    GROUP == "ControlSystemsBase" || dev_subpkg(GROUP) # This relies on the groups having the same name as the packages
    subpkg_path = joinpath(dirname(@__DIR__), "lib", GROUP)
    Pkg.test(PackageSpec(name = GROUP, path = subpkg_path))
end
