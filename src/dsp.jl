ControlSystems.tf(p::DSP.PolynomialRatio, h::Real = 1) = tf(coefb(p), coefa(p), h)

DSP.PolynomialRatio(G::TransferFunction) = DSP.PolynomialRatio(numpoly(G)[], denpoly(G)[])

function ControlSystems.TransferFunction(b::DSP.Biquad, h::Real = 1)
    b0, b1, b2, a1, a2 = b.b0, b.b1, b.b2, b.a1, b.a2
    tf([b0, b1, b2], [1, a1, a2], h)
end


"""
    Gs, k = seriesform(G::TransferFunction)

Convert a transfer function `G` to a vector of second-order transfer functions, the product of which equals `G`. `k` is a scalar gain.
"""
function seriesform(G::TransferFunction{<:Discrete})
    Gs = DSP.SecondOrderSections(DSP.PolynomialRatio(G))
    bqs = TransferFunction.(Gs.biquads, G.Ts)
    bqs, Gs.g
end


function DSP.SecondOrderSections(G::TransferFunction{<:Discrete})
    DSP.SecondOrderSections(DSP.PolynomialRatio(G))
end

function ControlSystems.zpk(p::DSP.ZeroPoleGain, h::Real)
    z,p,k = p.z, p.p, p.z
    zpk(z,p,k,h)
end

using Test
G = DemoSystems.resonant()*DemoSystems.resonant(ω0=2) |> tf
Gd = c2d(G, 0.1)
Gs, k = seriesform(Gd)
@test k*prod(Gs) ≈ Gd
Gds = DSP.SecondOrderSections(Gd)
