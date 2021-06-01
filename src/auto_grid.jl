

"""
values, grid = auto_grid(func, lims::Tuple, xscale, yscales; start_gridpoints, maxgridpoints, mineps, reltol_dd, abstol_dd)
    
    Arguments:
    Compute a `grid` that tries to capture all features of the function `func` in the range `lim=(xmin,xmax)`.
    `func`: `Function` so that that returns a `Tuple` of values, e.g. `func(x) = (v1,v2,...,vN)` where each `vi` is a scalar or AbstractArray.
        note that the behavior for non scalar vi might be unintuitive.
    The algorithm will try to produce a grid so that features are captured in all of the values, but widely varying scaled might cause problems.
    `xscale`: `Tuple` containing a monotone function and its inverse corresponding to the default axis.
        Example: If you want to produce a plot with logarithmic x-axis, then you should set `xscale=(log10, exp10)`.
    `yscales`: `Tuple` containing one monotone function for each of the values `(v1,v2,...,vN)` corresponding to the default axis.
    
    Kwargs:
    `start_gridpoints`: The number of initial grid points.
    `maxgridpoints`: A warning will be thrown if this number of gridpoints are reached.
    `mineps`: Lower limit on how small steps will be taken (in `xscale`).
        Example: with `xscale=(log10, exp10)` and `eps=1e-2`then `log10(grid[k+1])-log10(grid[k]) > (1e-2)/2`, i.e `grid[k+1] > 1.012grid[k]`.
    `reltol_dd = 0.05`, `abstol_dd = 0.01`: The criteria for further refinement is
        `norm(2*y2-y1-y3) < max(reltol_dd*max(norm(y2-y1), norm(y3-y2)), abstol_dd)`,
        where `y1,y2,y3` (in yscale) are 3 equally spaced points (in xscale).

    Output: 
    `values`: `Tuple` with one vector per output value of `func`.
    `grid`: `Vector` with the grid corresponding to the values, i.e `func(grid[i])[j] = values[j][i]`

    Example: The following code computes a grid and plots the functions: abs(sinc(x)) and cos(x)+1
        in lin-log and log-log plots respectively.
    
    func = x -> (abs(sinc(x)),cos(x)+1)
    lims = (0.1,7.0)
    xscale = (log10, exp10)
    yscales = (identity, log10)
    y, x = auto_grid(func, lims, xscale, yscales, mineps=1e-3)
    
    plot(x, y[1], m=:o; xscale=:log10, layout=(2,1), yscale=:identity, subplot=1, lab="sinc(x)", size=(800,600))
    plot!(x, y[2], m=:o, xscale=:log10, subplot=2, yscale=:log10, lab="cos(x)+1", ylims=(1e-5,10))

"""
function auto_grid(func, lims::Tuple, xscale::Tuple, yscales::Tuple;
                start_gridpoints=30, maxgridpoints=2000,
                mineps=(xscale[1](lims[2])-xscale[1](lims[1]))/10000,
                reltol_dd=0.05, abstol_dd=1e-2)
    
    # Linearly spaced grid in xscale
    init_grid = LinRange(xscale[1](lims[1]), xscale[1](lims[2]), start_gridpoints)
    # Current count of number mindpoints
    num_gridpoints = start_gridpoints

    # Current left point
    xleft = init_grid[1]
    xleft_identity = xscale[2](xleft)
    leftvalue = func(xleft_identity)
    # The full set of gridpoints
    grid = [xleft_identity,]
    # Tuple with list of all values (Faster than list of tuples + reshaping)
    values = ntuple(i -> [leftvalue[i],], length(yscales))
    
    for xright in init_grid[2:end]
        # Scale back to identity
        xright_identity = xscale[2](xright)
        rightvalue = func(xright_identity)
        # Refine (if nessesary) section (xleft,xright)
        num_gridpoints = refine_grid!(values, grid, func, xleft, xright, leftvalue, rightvalue, xscale, yscales, maxgridpoints, mineps, num_gridpoints; reltol_dd, abstol_dd)
        # We are now done with [xleft,xright]
        # Continue to next segment
        xleft = xright
        leftvalue = rightvalue
    end
    return values, grid
end

function push_multiple!(vecs::Tuple, vals::Tuple)
    push!.(vecs, vals)
end

# Given range (xleft, xright) potentially add more gridpoints. xleft is assumed to be already added, xright will be added
function refine_grid!(values, grid, func, xleft, xright, leftvalue, rightvalue, xscale, yscales, maxgridpoints, mineps, num_gridpoints; reltol_dd=0.05, abstol_dd=1e-2)
    # In scaled grid
    midpoint = (xleft+xright)/2
    # Scaled back
    midpoint_identity = xscale[2](midpoint)
    midvalue = func(midpoint_identity)

    (num_gridpoints >= maxgridpoints) && @warn "Maximum number of gridpoints reached in refine_grid! at $midpoint_identity, no further refinement will be made. Increase maxgridpoints to get better accuracy." maxlog=1
    #mineps in scaled version, abs to avoid assuming monotonly increasing scale
    if (abs(xright - xleft) >= mineps) && (num_gridpoints < maxgridpoints) && !is_almost_linear(xleft, midpoint, xright, leftvalue, midvalue, rightvalue, xscale, yscales; reltol_dd, abstol_dd)
        num_gridpoints = refine_grid!(values, grid, func, xleft,    midpoint, leftvalue, midvalue,   xscale, yscales, maxgridpoints, mineps, num_gridpoints; reltol_dd, abstol_dd)
        num_gridpoints = refine_grid!(values, grid, func, midpoint, xright,   midvalue,  rightvalue, xscale, yscales, maxgridpoints, mineps, num_gridpoints; reltol_dd, abstol_dd)
    else
        # No more refinement needed, add computed points
        push!(grid, midpoint_identity)
        push!(grid, xscale[2](xright))
        push_multiple!(values, midvalue)
        push_multiple!(values, rightvalue)
        num_gridpoints += 2
    end
    return num_gridpoints
end

# TODO We can scale when saving instead of recomputing scaling
# Svectors should also be faster in general
function is_almost_linear(xleft, midpoint, xright, leftvalue, midvalue, rightvalue, xscale, yscales; reltol_dd=0.05, abstol_dd=1e-2)
    # We assume that x2-x1 \approx  x3-x2, so we need to check that y2-y1 approx y3-y2

    # x1 = xscale.(xleft)
    # x2 = xscale.(midpoint)
    # x3 = xscale.(xright)
    y1 = apply_tuple2(yscales, leftvalue)
    y2 = apply_tuple2(yscales, midvalue)
    y3 = apply_tuple2(yscales, rightvalue)
    d1 = (y2-y1)
    d2 = (y3-y2)
    # Essentially low second derivative compared to derivative, so that linear approximation is good.
    # Second argument to avoid too small steps when derivatives are small
    norm(d1-d2) < max(reltol_dd*max(norm(d1), norm(d2)), abstol_dd)
end

# Seems to be fast enough
apply_tuple2(fs::Tuple, x) = [fs[i].(x[i]) for i in 1:length(fs)]