const Polynomials = ControlSystemsBase.Polynomials
import ControlSystemsBase.RootLocusResult
@userplot Rlocusplot


function getpoles(G, K::Number; kwargs...)
    issiso(G) || error("root locus only supports SISO systems")
    G isa TransferFunction || (G = tf(G))
    P = numpoly(G)[]
    Q = denpoly(G)[]
    T = float(typeof(K))
    ϵ = eps(T)
    nx = length(Q)
    D = zeros(nx-1, nx-1) # distance matrix
    prevpoles = ComplexF64[]
    temppoles = zeros(ComplexF64, nx-1)
    f = function (_,_,k)
        if k == 0 && length(P) > length(Q)
            # More zeros than poles, make sure the vector of roots is of correct length when k = 0
            # When this happens, there are fewer poles for k = 0, these poles can be seen as being located somewhere at Inf
            # We get around the problem by not allowing k = 0 for non-proper systems.
            k = ϵ
        end
        newpoles = ComplexF64.(Polynomials.roots(k[1]*P+Q))
        if !isempty(prevpoles)
            D .= abs.(newpoles .- transpose(prevpoles))
            assignment, cost = Hungarian.hungarian(D)
            for i = 1:nx-1
                temppoles[assignment[i]] = newpoles[i]
            end
            newpoles .= temppoles
        end
        prevpoles = newpoles
        newpoles
    end
    prob       = OrdinaryDiffEq.ODEProblem(f,f(0.,0.,0),(0,K[end]))
    integrator = OrdinaryDiffEq.init(prob,OrdinaryDiffEq.Tsit5(),reltol=1e-8,abstol=1e-8)
    ts         = Vector{Float64}()
    poleout    = Vector{Vector{ComplexF64}}()
    push!(poleout,integrator.k[1])
    push!(ts,0)
    for i in integrator
        push!(poleout,integrator.k[end])
        push!(ts,integrator.t[1])
    end
    poleout = copy(hcat(poleout...)')
    poleout, ts
end


function getpoles(G, K::AbstractVector{T}) where {T<:Number}
    issiso(G) || error("root locus only supports SISO systems")
    G isa TransferFunction || (G = tf(G))
    P, Q = numpoly(G)[], denpoly(G)[]
    poleout = Matrix{ComplexF64}(undef, Polynomials.degree(Q), length(K))
    nx = length(Q)
    D = zeros(nx-1, nx-1) # distance matrix
    temppoles = zeros(ComplexF64, nx-1)
    for (i, k) in enumerate(K)
        k == 0 && length(P) > length(Q) && (k = eps(T))
        poleout[:,i] = ComplexF64.(Polynomials.roots(k[1]*P+Q))
        if i > 1
            D .= abs.(poleout[:,i] .- transpose(poleout[:,i-1]))
            assignment, _ = Hungarian.hungarian(D)
            foreach(k->temppoles[assignment[k]] = poleout[:,i][k], 1:nx-1)
            poleout[:,i] .= temppoles
        end
    end
    copy(poleout'), K
end

"""
    getpoles(sys::StateSpace, K::AbstractMatrix; tol = 1e-2, initial_stepsize = 1e-3, output=false)

Compute the poles of the closed-loop system defined by `sys` with feedback gains `γ*K` where `γ` is a scalar that ranges from 0 to 1.

If `output = true`, `K` is assumed to be an output feedback matrix of dim `(nu, ny)`
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

    while k_scalar < 1.0
        # Propose a new k_scalar value
        next_k_scalar = min(1.0, k_scalar + stepsize)

        # Calculate poles for the proposed next_k_scalar
        A_cl_proposed = A - next_k_scalar * B * K_matrix
        current_poles_proposed = eigvals(A_cl_proposed)

        # Calculate cost using Hungarian algorithm
        D = zeros(nx, nx)
        for r in 1:nx
            for c in 1:nx
                D[r, c] = abs(current_poles_proposed[r] - prevpoles[c])
            end
        end
        assignment, cost = Hungarian.hungarian(D)

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

    return hcat(poleout_list...)' |> copy, k_scalars_collected .* Ref(K_matrix) # Return transposed pole matrix and k_values
end


"""
    roots, Z, K = rlocus(P::LTISystem, K = 500)

Compute the root locus of the SISO LTISystem `P` with a negative feedback loop and feedback gains between 0 and `K`. `rlocus` will use an adaptive step-size algorithm to determine the values of the feedback gains used to generate the plot.

`roots` is a complex matrix containing the poles trajectories of the closed-loop `1+k⋅G(s)` as a function of `k`, `Z` contains the zeros of the open-loop system `G(s)` and `K` the values of the feedback gain.

If `K` is a matrix and `P` a `StateSpace` system, the poles are computed as `K` ranges from `0*K` to `1*K`. In this case, `K` is assumed to be a state-feedback matrix of dimension `(nu, nx)`. To compute the poles for output feedback, use, pass `output = true` and `K` of dimension `(nu, ny)`.
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
        redata[1,:], imdata[1,:]
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

Plot the root locus of the SISO LTISystem `P` as computed by `rlocus`.
"""
@recipe function rlocusplot(::Type{Rlocusplot}, p::Rlocusplot; K=500, output=false)
    if length(p.args) >= 2
        rlocus(p.args[1], p.args[2]; output)
    else
        rlocus(p.args[1]; K=K)
    end
end