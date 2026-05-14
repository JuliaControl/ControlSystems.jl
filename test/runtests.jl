using Pkg
using Test, LinearAlgebra, Random
using Aqua

# Helper to call dev on one of the packages in folder /lib
function dev_subpkg(subpkg)
    subpkg_path = joinpath(dirname(@__DIR__), "lib", subpkg)
    Pkg.develop(PackageSpec(path = subpkg_path))
end

dev_subpkg("ControlSystemsBase") # Always dev this package to test with the latest code

const GROUP = get(ENV, "GROUP", "All") # Get the GROUP attribute from the test.yml file, default to "All" for testing locally

@show GROUP

if GROUP ∈ ("ControlSystems", "All")
    include("runtests_controlsystems.jl")
end

if GROUP == "All"
    for GROUP in readdir(joinpath(dirname(@__DIR__), "lib")) # Loop over all subpackages
        subpkg_path = joinpath(dirname(@__DIR__), "lib", GROUP)
        Pkg.test(PackageSpec(name = GROUP, path = subpkg_path), coverage = true)
    end
    @testset "Aqua" begin
        import ControlSystems
        Aqua.test_all(ControlSystems;
            ambiguities = false, # causes 100s of hits in all dependencies
            stale_deps = true,
        )
    end

elseif GROUP != "ControlSystems"
    # dev_subpkg(GROUP) # Do this if more sub packages are added, don't forget to avoid doing it if GROUP is CSBase
    subpkg_path = joinpath(dirname(@__DIR__), "lib", GROUP)
    Pkg.test(PackageSpec(name = GROUP, path = subpkg_path), coverage = true)
end

