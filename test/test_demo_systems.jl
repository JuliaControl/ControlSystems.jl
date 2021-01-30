@testset "test_demosystems" begin

# Just some very simple tests of systemdepot
@test size(DemoSystems.woodberry()) == (2,2)
@test size(DemoSystems.fotd()) == (1,1)
@test freqresp(DemoSystems.fotd(τ=2, T=5), [0.0])[:] == [1.0]


Random.seed!(10)

sys = ssrand(1,proper=true)
@test size(sys) == (1,1)
@test ControlSystems.nstates(sys) == 2
@test iszero(sys.D)

@test all(isstable, ssrand(1,1,2,stable=true) for _=1:100)
@test all(isstable, ssrand(1,1,2,stable=true,Ts=0.1) for _=1:100)
@test all(isstable, ssrand(1,1,2,stable=true,Ts=10) for _=1:100)

sys = ssrand(2,2,5,proper=false,stable=true)
@test size(sys) == (2,2)
@test ControlSystems.nstates(sys) == 5
@test isstable(sys)
@test !iszero(sys.D)

Random.seed!(20)
sys = ssrand(2,2,5,stable=false)
@test !isstable(sys)

sys = ssrand(2,2,5,proper=false,stable=true, Ts=0.5)
@test size(sys) == (2,2)
@test ControlSystems.nstates(sys) == 5
@test isstable(sys)
@test !iszero(sys.D)

sys = ssrand(2,2,5,proper=false,stable=false, Ts=0.5)
@test !isstable(sys)



c = [2.5]
sysd = DemoSystems.ssfir(c)
y = impulse(sysd, 10)[1]
@test y[1] == c[1]
@test y[2:11] == zeros(10)


c = [0.1, 2.5]
sysd = DemoSystems.ssfir(c)
y = impulse(sysd, 10)[1]
@test y[1:2] == c[1:2]
@test y[3:11] == zeros(9)


c = [1,2,5,3]
sysd = DemoSystems.ssfir(c)
y = impulse(sysd, 10)[1]
@test y[1:4] == c
@test y[5:11] == zeros(7)

end
