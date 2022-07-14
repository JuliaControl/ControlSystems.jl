"""
    ratelimit(val; Tf)
    ratelimit(lower, upper; Tf)

Create a nonlinearity that limits the rate of change of a signal, roughly equivalent to `` 1/s ∘ sat ∘ s``. `Tf` controls the filter time constant on the derivative used to calculate the rate. 
$nonlinear_warning
"""
function ratelimit(args...; Tf)
    T = promote_type(typeof(Tf), typeof.(args)...)
    # Filtered derivative
    A = [0 1; -2/Tf^2 -2/Tf]
    B = [0; 1]
    C = [4/Tf^2 0]
    F = ss(A,B,C,0)

    integrator = tf([T(Tf), one(T)], [one(T), zero(T)])
    integrator*saturation(args...)*F
end

