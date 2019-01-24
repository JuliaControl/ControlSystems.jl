@testset "test_arx" begin
N = 20
t = 1:N
u = randn(N)
G = tf(0.8, [1,-0.9], 1)
y = lsim(G,u,t)[1][:]

na,nb = 1,1
yr,A = getARXregressor(y,u,na,nb)
@test length(yr) == N-na
@test size(A) == (N-na, na+nb)

@test yr == y[na+1:end]
@test A[:,1] == y[1:end-na]
@test A[:,2] == u[1:end-1]

na = 2
Gh,Σ = arx(1,y,u,na,nb)
@test Gh ≈ G # Should recover the original transfer function exactly
ω=exp10.(range(-2, stop=1, length=200))
# bodeplot(G,ω)
# bodeconfidence!(Gh, Σ, ω=ω, color=:blue, lab="Est")

# Test MISO estimation
u2 = randn(N)
G2 = [G tf(0.5, [1, -0.9], 1)]
y2 = lsim(G2,[u u2],t)[1][:]

nb = [1,1]
Gh2,Σ = arx(1,y2,[u u2],na,nb)

@test Gh2 ≈ G2

end
