const Polynomials = ControlSystemsBase.Polynomials
@userplot Rlocusplot
@deprecate rlocus(args...;kwargs...) rlocusplot(args...;kwargs...)


function getpoles(G, K)
    issiso(G) || error("root locus only supports SISO systems")
    G isa TransferFunction || (G = tf(G))
    P = numpoly(G)[]
    Q = denpoly(G)[]
    T = float(eltype(K))
    ϵ = eps(T)
    nx = length(Q)
    D = zeros(nx-1, nx-1) # distance matrix
    prevpoles = ComplexF64[]
    temppoles = zeros(ComplexF64, nx-1)
    f = function (y,_,k)
        if k == 0 && length(P) > length(Q) 
            # More zeros than poles, make sure the vector of roots is of correct length when k = 0
            # When this happens, there are fewer poles for k = 0, these poles can be seen as beeing located somewhere at Inf
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
    for i in integrator
        push!(poleout,integrator.k[1])
        push!(ts,integrator.t[1])
    end
    poleout = hcat(poleout...)'
    poleout, ts
end



"""
    rlocusplot(P::LTISystem; K)

Computes and plots the root locus of the SISO LTISystem P with
a negative feedback loop and feedback gains between 0 and `K`. `rlocusplot` will use an adaptive step-size algorithm to
determine the values of the feedback gains used to generate the plot.
"""
rlocusplot
@recipe function rlocusplot(p::Rlocusplot; K=500)
    P = p.args[1]
    K = K isa Number ? range(1e-6,stop=K,length=10000) : K
    Z = tzeros(P)
    roots, K = getpoles(P,K)
    redata = real.(roots)
    imdata = imag.(roots)
    
    ylims --> (max(-50,minimum(imdata) - 1), min(50,maximum(imdata) + 1))
    xlims --> (max(-50,minimum(redata) - 1), clamp(maximum(redata) + 1, 1, 50))
    framestyle --> :zerolines
    title --> "Root locus"
    xguide --> "Re(roots)"
    yguide --> "Im(roots)"
    form(k, p) = Printf.@sprintf("%.4f", k) * "  pole=" * Printf.@sprintf("%.3f%+.3fim", real(p), imag(p))
    @series begin
        legend --> false
        hover := "K=" .* form.(K,roots)
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
end