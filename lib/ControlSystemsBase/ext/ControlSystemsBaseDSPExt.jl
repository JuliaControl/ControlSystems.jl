module ControlSystemsBaseDSPExt
using ControlSystemsBase
using ControlSystemsBase: issiso, numvec, denvec
import ControlSystemsBase: TransferFunction, seriesform, zpk, tf
import DSP

tf(p::DSP.PolynomialRatio{:z}, h::Real = 1) = tf(DSP.coefb(p), DSP.coefa(p), h)
tf(p::DSP.PolynomialRatio{:s}) = tf(DSP.coefb(p), DSP.coefa(p))
tf(p::DSP.ZeroPoleGain{:z}, h::Real = 1) = tf(DSP.PolynomialRatio(p), h)
tf(p::DSP.ZeroPoleGain{:s}) = tf(DSP.PolynomialRatio(p))

function DSP.PolynomialRatio(G::TransferFunction{<:Discrete})
    DSP.PolynomialRatio(
        numvec(G)[],
        denvec(G)[],
    )
end

function TransferFunction(b::DSP.Biquad, h::Real = 1)
    b0, b1, b2, a1, a2 = b.b0, b.b1, b.b2, b.a1, b.a2
    tf([b0, b1, b2], [1, a1, a2], h)
end


function seriesform(G::TransferFunction{<:Discrete})
    Gs = DSP.SecondOrderSections(DSP.PolynomialRatio(G))
    bqs = TransferFunction.(Gs.biquads, G.Ts)
    bqs, Gs.g
end


function DSP.SecondOrderSections(G::TransferFunction{<:Discrete})
    DSP.SecondOrderSections(DSP.PolynomialRatio(G))
end

function zpk(p::DSP.ZeroPoleGain, h::Real)
    z,p,k = p.z, p.p, p.k
    zpk(z, p, k, h)
end

"""
    DSP.filt(P::ControlSystemsBase.TransferFunction, u)

Use a transfer function `P` to filter a signal `u`. This is equivalent to `lsim(P, u').y[:]`, but may be more efficient for single-input, single-output systems when the state sequence is not needed.
"""
function DSP.filt(P::ControlSystemsBase.TransferFunction, u, args...)
    issiso(P) || error("Only single-input, single-output systems are supported in filt, call lsim instead.")
    b, a = numvec(P)[], denvec(P)[]
    nb, na = length(b), length(a)
    if nb <= na
        b = [zeros(na - nb); b]
    end
    DSP.filt(b, a, u, args...)
end


function DSP.filtfilt(P::ControlSystemsBase.TransferFunction, u, args...)
    issiso(P) || error("Only single-input, single-output systems are supported in filtfilt.")
    b, a = numvec(P)[], denvec(P)[]
    nb, na = length(b), length(a)
    if nb <= na
        b = [zeros(na - nb); b]
    end
    DSP.filtfilt(b, a, u, args...)
end


end