# MonteCarloMeasurements depends on GenericSchur, so this file must run after the tests that require generic linear algebra to be unavailable.
import MonteCarloMeasurements as MCM
@testset "MonteCarloMeasurements interoperability" begin
    @info "Testing MonteCarloMeasurements interoperability"
    # The roots of a polynomial with particle coefficients must be ordered consistently
    # over the particles so that the particle-valued roots form conjugate pairs.
    ω = MCM.:..(4, 12)
    ζ = MCM.:..(0.054, 0.084)
    G = tf([2*ζ/ω, 1], [1/ω^2, 2*ζ/ω, 1, 0, 0])
    C = tf([1.0, 1], [0.1, 1])
    Gc = feedback(C*G, 1)
    Z = zpk(Gc)
    p = Z.matrix[1].p
    @test length(p) == 5
    @test count(x -> !iszero(imag(x)), p) == 4
    w = exp10.(range(-1, 2, length=20))
    @test MCM.pmean.(bode(Z, w)[1]) ≈ MCM.pmean.(bode(Gc, w)[1]) rtol=1e-6
end
