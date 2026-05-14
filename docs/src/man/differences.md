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
- Julia compiles code just before it is called the first time. This introduces a noticeable lag the first time a function is called, and can make packages take a while to load.


If you find other noteworthy differences between ControlSystems.jl and other pieces of control-related software, please consider [opening an issue on the github repository](https://github.com/JuliaControl/ControlSystems.jl/issues/new?template=Blank+issue).