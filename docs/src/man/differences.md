# Noteworthy Differences from other Languages
If you are new to the Julia programming language, you are encouraged to visit the documentation page on [noteworthy differences between Julia and other programming languages](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).

The rest of this page will list noteworthy differences between ControlSystems.jl and other pieces of control-systems software.

- Functions to calculate poles and zeros of systems are named using their plural forms, i.e., [`poles`](@ref) instead of `pole`, and [`tzeros`](@ref) instead of `tzero`.
- Simulation using [`lsim`](@ref), [`step`](@ref), [`impulse`](@ref) returns arrays where time is in the second dimension rather than in the first dimension (applies also to `freqresp, bode, nyquist` etc.). Julia uses a *column major* memory layout, and this choice is made for performance reasons.
- Functions [`are`](@ref), [`lqr`](@ref) and [`kalman`](@ref) have slightly different signatures in julia compared to in other languages. More advanced LQG functionalities are located in [RobustAndOptimalControl.jl](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/).
- Simulation using [`lsim`](@ref), [`step`](@ref), [`impulse`](@ref) etc. return a structure that can be plotted. These functions never plot anything themselves.
- Functions [`bode`](@ref), [`nyquist`](@ref) etc. never produce a plot. Instead, see [`bodeplot`](@ref), [`nyquistplot`](@ref) etc.
- In Julia, functionality is often split up into several different packages. You may therefore have to install and use additional packages in order to cover all your needs. See [Ecosystem](@ref) for a collection of control-related packages.
- In Julia, `1` has a different type than `1.0`, and the types in ControlSystemsBase.jl respect the types chosen by the user. As an example, `tf(1, [1, 1])` is a transfer function with integer coefficients, while `tf(1.0, [1, 1])` will promote all coefficients to `Float64`.
- Julia treats matrices and vectors as different types, in particular, column vectors and row vectors are not interchangeable. 
- In Julia, code can often be differentiated using automatic differentiation. When using ControlSystems.jl, we recommend trying [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl/) for AD. An example making use of this is available in [Automatic Differentiation](@ref).
- In Julia, the source code is often very readable. If you want to learn how a function is implemented, you may find the macros `@edit` or `@less` useful to locate the source code.
- If you run into an issue (bug) with a Julia package, you can share this issue (bug report) on the package's github page and it will often be fixed promptly. To open an issue with ControlSystems.jl, [click here](https://github.com/JuliaControl/ControlSystems.jl/issues/new/choose). Thank you for helping out improving open-source software!
- Julia compiles code just before it is called the first time. This introduces a noticeable lag, and can make packages take a long time to load. If you want to speed up the loading of ControlSystems.jl, consider building a system image that includes ControlSystems.jl using [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/stable/). More info about this is available below under [Precompilation for faster load times](@ref)


If you find other noteworthy differences between ControlSystems.jl and other pieces of control-related software, please consider submitting a pull request (PR) to add to the list above. You can submit a PR by clicking on "Edit on GitHub" at the top of this page and then clicking on the icon that looks like a pen above the file viewer. A two-minute video on this process is available below
```@raw html
<iframe style="height: 315px; width: 560px" src="https://www.youtube.com/embed/ZpH1ry8qqfw" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
```

## Precompilation for faster load times
In order to make it faster to load the ControlSystems.jl package, you may make use of [PackageCompiler.jl](https://julialang.github.io/PackageCompiler.jl/stable/). 

!!! warning "For developers"
    If you intend to develop ControlSystem.jl, i.e., modify the source code, it's not recommended to build the package into the system image. We then recommend to build OrdinaryDiffEq into the system image since this package contributes the largest part of the loading time.

Building a custom system image can reduce the time to get started in a new Julia session, as an example:

- Without system image:
```julia
julia> @time using ControlSystems
  1.646961 seconds (2.70 M allocations: 173.558 MiB, 1.08% gc time, 2.06% compilation time)
```

- With OrdinaryDiffEq and Plots in the system image:
```julia
julia> @time using ControlSystems
  0.120975 seconds (413.37 k allocations: 27.672 MiB, 1.66% compilation time)
```


To build a system image with ControlSystems, save the following script in a file, e.g., `precompile_controlsystems.jl` (feel free to add any additional packages you may want to load).
```julia
using OrdinaryDiffEq # Include this if you want to develop ControlSystems.jl
using ControlSystems # Include this if you only want to use ControlSystems.jl
using Plots # In case you also want to use plotting functions

# Run some statements to make sure these are precompiled. Do not include this if you want to develop ControlSystems.jl
for P = StateSpace[ssrand(2,2,2), ssrand(2,2,2, Ts=0.1)]
    bodeplot(P)
    nyquistplot(P)
    plot(step(P, 10))
end
```

Then run the following
```julia
using PackageCompiler
PackageCompiler.create_sysimage(
    [
        :OrdinaryDiffEq,
        :Plots,
        :ControlSystems,
    ];
    precompile_execution_file = "precompile_execution_file",
    sysimage_path = "sys_ControlSystems_$(VERSION).so",
)
exit()
```

When you have created a system image, start Julia with the -J flag pointing to the system image that was created, named `sys_ControlSystems_<VERSION>.so`, [more details here](https://julialang.github.io/PackageCompiler.jl/dev/sysimages.html#Creating-a-sysimage-using-PackageCompiler). After this, loading the package should be very fast.

!!! warning "Updating packages"
    When you update installed julia packages, the update will not be reflected in the system image until the image is rebuilt. 

You can make vscode load this system image as well by adding
```json
"julia.additionalArgs": [
    "-J/path_to_sysimage/sys_ControlSystems_<VERSION>.so"
],
```
to `settings.json`.