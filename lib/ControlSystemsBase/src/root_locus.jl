import ControlSystemsBase.RootLocusResult
@userplot Rlocusplot

"""
    rlocusplot(siso_sys)
    rlocusplot(sys::StateSpace, K::Matrix; output=false)

Plot the root locus of a system under feedback.

If a SISO system is passed, the feedback gain `K` is a scalar that ranges from 0 to `K` (if provided). If a `StateSpace` system is passed, `K` is a matrix that defines the feedback gain, and the poles are computed as `K` ranges from `0*K` to `1*K`. In this case, `K` is assumed to be a state-feedback matrix of dimension `(nu, nx)`. To compute the poles for output feedback, pass `output = true` and `K` of dimension `(nu, ny)`.
"""
rlocusplot

"""
    pole_reference_scale(poles...)

Return a magnitude below which the step-size control of the root locus measures the displacement of a pole in absolute rather than relative terms. The scale is a small fraction of the smallest nonzero magnitude among the poles and zeros passed in, where magnitudes below `sqrt(eps())` times the largest magnitude are regarded as zero since they are typically the result of rounding errors in the computation of a root at the origin. This bounds the number of steps required for a pole that approaches or departs from the origin.
"""
function pole_reference_scale(poles...)
    magnitudes = [abs(x) for p in poles for x in p if isfinite(x)]
    largest = maximum(magnitudes; init=0.0)
    smallest = minimum(a for a in magnitudes if a > sqrt(eps())*largest; init=Inf)
    isfinite(smallest) ? 1e-3*smallest : 1.0
end

"""
    relative_step_cost(D, assignment, prevpoles, s_ref)

Return the sum over all poles of the distance moved in one step, each distance divided by the magnitude of the pole before the step (but at least `s_ref`). Measuring the displacement relative to the pole magnitude makes the number of steps grow logarithmically rather than linearly with the distance the poles travel.
"""
function relative_step_cost(D, assignment, prevpoles, s_ref)
    cost = 0.0
    for (r, c) in enumerate(assignment)
        cost += D[r, c] / max(abs(prevpoles[c]), s_ref)
    end
    cost
end

"""
    improper_initial_gain(P, Q, K; radius_factor = 10, rtol = 1e-2, npoints = 64)

Return the gain at which the root locus of the improper transfer function `P/Q` is started, at most `K/10`. At this gain, the `degree(P) - degree(Q)` poles that originate at infinity have a magnitude larger than `radius_factor` times the largest magnitude among the open-loop poles and zeros, and each of the remaining poles is located within a distance `rtol*|p|` of an open-loop pole `p` (or a cluster thereof).

The bound follows from Rouché's theorem: if `k|P(s)| < |Q(s)|` on a closed contour, `Q + kP` has as many roots inside the contour as `Q`. The contour condition is evaluated at `npoints` points on each contour.
"""
function improper_initial_gain(P, Q, K; radius_factor = 10, rtol = 1e-2, npoints = 64)
    ol_poles = Polynomials.roots(Q)
    largest = maximum(abs, [ol_poles; Polynomials.roots(P)]; init=0.0)
    largest > 0 || (largest = 1.0) # All poles and zeros are located at the origin, the locus has no intrinsic scale
    circle = cis.(range(0, 2π, length=npoints+1)[1:end-1])
    contour = radius_factor*largest .* circle
    for p in ol_poles
        append!(contour, p .+ rtol*max(abs(p), rtol*largest) .* circle)
    end
    k0 = minimum(s->abs(Q(s))/abs(P(s)), contour)
    k0 = min(k0, K/10)
    k0 > 0 ? k0 : eps(float(K))
end


function getpoles(G, K::Number; tol = 1e-2, initial_stepsize = 1e-3, kwargs...)
    issiso(G) || error("root locus with scalar gain only supports SISO systems, did you intend to pass a feedback gain matrix `K`?")
    G isa TransferFunction || (G = tf(G))
    P = numpoly(G)[]
    Q = denpoly(G)[]
    T = float(typeof(K))
    ϵ = eps(T)
    improper = Polynomials.degree(P) > Polynomials.degree(Q)
    npoles = max(Polynomials.degree(P), Polynomials.degree(Q)) # Number of closed-loop poles for k > 0
    
    # Scale tolerance with system order
    tol = tol * npoles
    
    poleout_list = Vector{Vector{ComplexF64}}() # To store pole sets at each accepted step
    k_scalars_collected = Float64[] # To store accepted k_scalar values
    
    prevpoles = ComplexF64[] # Initialize prevpoles for the first iteration
    temppoles = zeros(ComplexF64, npoles)
    D = zeros(npoles, npoles) # distance matrix
    
    stepsize = initial_stepsize
    # For an improper system, the closed-loop system has more poles than the open-loop system for k > 0, the additional poles originate at infinity. The locus is thus started at a gain k > 0 at which these poles have a large magnitude.
    k_scalar = improper && K > 0 ? improper_initial_gain(P, Q, K) : 0.0
    
    # Function to compute poles for a given k value
    compute_poles = function(k)
        if k == 0 && improper
            # Make sure the vector of roots is of correct length, the additional poles are located at a large magnitude
            k = ϵ
        end
        ComplexF64.(Polynomials.roots(k*P+Q))
    end
    
    initial_poles = compute_poles(k_scalar)
    push!(poleout_list, initial_poles)
    push!(k_scalars_collected, k_scalar)
    prevpoles = initial_poles # Set prevpoles for the first actual step
    s_ref = pole_reference_scale(Polynomials.roots(Q), compute_poles(K), Polynomials.roots(P))
    
    while k_scalar < K
        # Propose a new k_scalar value
        next_k_scalar = min(K, k_scalar + stepsize)
        
        # Calculate poles for the proposed next_k_scalar
        current_poles_proposed = compute_poles(next_k_scalar)
        
        # Calculate cost using Hungarian algorithm
        if !isempty(prevpoles)
            D .= abs.(current_poles_proposed .- transpose(prevpoles))
            assignment, _ = Hungarian.hungarian(D)
            cost = relative_step_cost(D, assignment, prevpoles, s_ref)
        else
            cost = 0.0
            assignment = collect(1:npoles)
        end
        
        # Adaptive step size logic
        if cost > 2 * tol # Cost is too high, reject step and reduce stepsize
            stepsize /= 2.0
            # Ensure stepsize doesn't become too small
            if stepsize < 100*eps(T)
                @warn "Step size became extremely small, potentially stuck. Breaking loop."
                break
            end
            # Do not update k_scalar, try again with smaller stepsize
        else # Step is acceptable
            # Sort poles using the assignment from Hungarian algorithm
            if !isempty(prevpoles)
                for i = 1:npoles
                    temppoles[assignment[i]] = current_poles_proposed[i]
                end
                current_poles_sorted = copy(temppoles)
            else
                current_poles_sorted = current_poles_proposed
            end
            
            # Accept the step
            push!(poleout_list, current_poles_sorted)
            push!(k_scalars_collected, next_k_scalar)
            prevpoles = current_poles_sorted # Update prevpoles for the next iteration
            k_scalar = next_k_scalar # Advance k_scalar
            
            if cost < tol # Cost is too low, increase stepsize
                stepsize *= 1.1
            end
            # Cap stepsize to prevent overshooting K significantly in a single step
            stepsize = min(stepsize, K/10)
        end
        
        # Break if k_scalar has reached or exceeded K
        if k_scalar >= K
            break
        end
    end
    
    return copy(transpose(reduce(hcat, poleout_list))), k_scalars_collected
end


function getpoles(G, K::AbstractVector{T}) where {T<:Number}
    issiso(G) || error("root locus with scalar gain only supports SISO systems, did you intend to pass a feedback gain matrix `K`?")
    G isa TransferFunction || (G = tf(G))
    P, Q = numpoly(G)[], denpoly(G)[]
    npoles = max(Polynomials.degree(P), Polynomials.degree(Q)) # Number of closed-loop poles for k > 0
    poleout = Matrix{ComplexF64}(undef, npoles, length(K))
    D = zeros(npoles, npoles) # distance matrix
    temppoles = zeros(ComplexF64, npoles)
    for (i, k) in enumerate(K)
        # For an improper system, the additional poles are located at infinity for k = 0, they are approximated by poles of large magnitude
        k == 0 && Polynomials.degree(P) > Polynomials.degree(Q) && (k = eps(float(T)))
        poleout[:,i] = ComplexF64.(Polynomials.roots(k[1]*P+Q))
        if i > 1
            D .= abs.(poleout[:,i] .- transpose(poleout[:,i-1]))
            assignment, _ = Hungarian.hungarian(D)
            foreach(k->temppoles[assignment[k]] = poleout[:,i][k], 1:npoles)
            poleout[:,i] .= temppoles
        end
    end
    copy(transpose(poleout)), K
end

"""
    getpoles(sys::StateSpace, K::AbstractMatrix; tol = 1e-2, initial_stepsize = 1e-3, output=false)

Compute the poles of the closed-loop system defined by `sys` with feedback gains `γ*K` where `γ` is a scalar that ranges from 0 to 1.

If `output = true`, `K` is assumed to be an output feedback matrix of dim `(nu, ny)`

The step size in `γ` is adapted such that the sum over all poles of the distance moved in one step, each relative to the magnitude of the pole, is at most `2tol` times the state dimension.
"""
function getpoles(sys::StateSpace, K_matrix::AbstractMatrix; tol = 1e-2, initial_stepsize = 1e-3, output=false)
    (; A, B, C) = sys
    nx = size(A, 1) # State dimension
    ny = size(C, 1) # Output dimension
    tol = tol*nx # Scale tolerance with state dimension
    # Check for compatibility of K_matrix dimensions with B
    if size(K_matrix, 2) != (output ? ny : nx)
        error("The number of columns in K_matrix ($(size(K_matrix, 2))) must match the state dimension ($(nx)) or output dimension ($(ny)) depending on whether output feedback is used.")
    end
    if size(K_matrix, 1) != size(B, 2)
        error("The number of rows in K_matrix ($(size(K_matrix, 1))) must match the number of inputs (columns of B, which is $(size(B, 2))).")
    end

    if output
        # We bake C into K here to avoid repeated multiplications below
        K_matrix = K_matrix * C
    end

    poleout_list = Vector{Vector{ComplexF64}}() # To store pole sets at each accepted step
    k_scalars_collected = Float64[] # To store accepted k_scalar values

    prevpoles = ComplexF64[] # Initialize prevpoles for the first iteration

    stepsize = initial_stepsize
    k_scalar = 0.0

    # Initial poles at k_scalar = 0.0
    A_cl_initial = A - 0.0 * B * K_matrix
    initial_poles = eigvals(A_cl_initial)
    push!(poleout_list, initial_poles)
    push!(k_scalars_collected, 0.0)
    prevpoles = initial_poles # Set prevpoles for the first actual step
    s_ref = pole_reference_scale(initial_poles, eigvals(A - B * K_matrix))
    D = zeros(nx, nx) # distance matrix

    while k_scalar < 1.0
        # Propose a new k_scalar value
        next_k_scalar = min(1.0, k_scalar + stepsize)

        # Calculate poles for the proposed next_k_scalar
        A_cl_proposed = A - next_k_scalar * B * K_matrix
        current_poles_proposed = eigvals(A_cl_proposed)

        # Calculate cost using Hungarian algorithm
        D .= abs.(current_poles_proposed .- transpose(prevpoles))
        assignment, _ = Hungarian.hungarian(D)
        cost = relative_step_cost(D, assignment, prevpoles, s_ref)

        # Adaptive step size logic
        if cost > 2 * tol # Cost is too high, reject step and reduce stepsize
            stepsize /= 2.0
            # Ensure stepsize doesn't become too small
            if stepsize < 100eps()
                @warn "Step size became extremely small, potentially stuck. Breaking loop."
                break
            end
            # Do not update k_scalar, try again with smaller stepsize
        else # Step is acceptable or too small
            # Sort poles using the assignment from Hungarian algorithm
            temppoles = zeros(ComplexF64, nx)
            for j = 1:nx
                temppoles[assignment[j]] = current_poles_proposed[j]
            end
            current_poles_sorted = temppoles

            # Accept the step
            push!(poleout_list, current_poles_sorted)
            push!(k_scalars_collected, next_k_scalar)
            prevpoles = current_poles_sorted # Update prevpoles for the next iteration
            k_scalar = next_k_scalar # Advance k_scalar

            if cost < tol # Cost is too low, increase stepsize
                stepsize *= 1.1
            end
            # Cap stepsize to prevent overshooting 1.0 significantly in a single step
            stepsize = min(stepsize, 1e-1)
        end

        # Break if k_scalar has reached or exceeded 1.0
        if k_scalar >= 1.0
            break
        end
    end

    return copy(transpose(reduce(hcat, poleout_list))), k_scalars_collected .* Ref(K_matrix) # Return transposed pole matrix and k_values
end


"""
    roots, Z, K = rlocus(P::LTISystem, K = 500)

Compute the root locus of the LTISystem `P` with a negative feedback loop and feedback gains between 0 and `K`. `rlocus` will use an adaptive step-size algorithm to determine the values of the feedback gains used to generate the plot.

`roots` is a complex matrix containing the poles trajectories of the closed-loop `1+k⋅G(s)` as a function of `k`, `Z` contains the zeros of the open-loop system `G(s)` and `K` the values of the feedback gain.

If `P` is improper, the closed-loop system has more poles than the open-loop system for `k > 0`, the additional poles originate at infinity. The root locus is then started at a gain `K[1] > 0` at which these poles have a large magnitude in relation to the open-loop poles and zeros.

If `K` is a matrix and `P` a `StateSpace` system, the poles are computed as `K` ranges from `0*K` to `1*K`. In this case, `K` is assumed to be a state-feedback matrix of dimension `(nu, nx)`. To compute the poles for output feedback, use, pass `output = true` and `K` of dimension `(nu, ny)`.

The keyword arguments `tol = 1e-2` and `initial_stepsize = 1e-3` control the adaptive step size. A step is accepted if the sum over all poles of the distance moved, each relative to the magnitude of the pole, is at most `2tol` times the number of poles. The number of steps thus grows logarithmically with the distance the poles travel, independent of the time unit of the system.
"""
function rlocus(P, K; kwargs...)
    Z = tzeros(P)
    roots, K = getpoles(P, K; kwargs...)
    ControlSystemsBase.RootLocusResult(roots, Z, K, P)
end

rlocus(P; K=500) = rlocus(P, K)


# This will be called on plot(rlocus(sys, args...))
@recipe function rootlocusresultplot(r::RootLocusResult)
    roots, Z, K = r
    array_K = eltype(K) <: AbstractArray
    redata = real.(roots)
    imdata = imag.(roots)
    all_redata = [vec(redata); real.(Z)]
    all_imdata = [vec(imdata); imag.(Z)]

    ylims --> (max(-50,minimum(all_imdata) - 1), min(50,maximum(all_imdata) + 1))
    xlims --> (max(-50,minimum(all_redata) - 1), clamp(maximum(all_redata) + 1, 1, 50))
    framestyle --> :zerolines
    title --> "Root locus"
    xguide --> "Re(roots)"
    yguide --> "Im(roots)"
    form(k, p) = Printf.@sprintf("%.4f", k) * "  pole=" * Printf.@sprintf("%.3f%+.3fim", real(p), imag(p))
    @series begin
        legend --> false
        if !array_K
            hover := "K=" .* form.(K,roots)
        end
        label := ""
        redata, imdata
    end
    @series begin
        seriestype := :scatter
        markershape --> :circle
        markersize --> 10
        label --> "Zeros"
        real.(Z), imag.(Z)
    end
    @series begin
        seriestype := :scatter
        markershape --> :xcross
        markersize --> 10
        label --> "Open-loop poles"
        ol_poles = poles(r.sys)
        real.(ol_poles), imag.(ol_poles)
    end
    if array_K
        @series begin
            seriestype := :scatter
            markershape --> :diamond
            markersize --> 10
            label --> "Closed-loop poles"
            redata[end,:], imdata[end,:]
        end
    end
end


"""
    rlocusplot(P::LTISystem; K)
    rlocusplot(P::StateSpace, K::Matrix; output = false)

Plot the root locus of the LTISystem `P` as computed by `rlocus`.
"""
@recipe function rlocusplot(::Type{Rlocusplot}, p::Rlocusplot; K=500, output=false)
    if length(p.args) >= 2
        rlocus(p.args[1], p.args[2]; output)
    else
        rlocus(p.args[1]; K=K)
    end
end
