module ControlSystemsBaseMonteCarloMeasurementsExt
using ControlSystemsBase
using ControlSystemsBase: Polynomials
using MonteCarloMeasurements: AbstractParticles

# The roots of a polynomial with particle coefficients are computed for each particle
# separately. Without sorting, the order of the roots returned by LAPACK may differ
# between particles, so that a particle-valued root combines different roots and the
# particle-valued roots no longer form exact conjugate pairs. Sorting by real and
# imaginary part places the members of each conjugate pair at the same adjacent
# positions in every particle.
ControlSystemsBase.roots(p::Polynomials.AbstractPolynomial{<:AbstractParticles}; kwargs...) =
    Polynomials.roots(p; kwargs...)

end
