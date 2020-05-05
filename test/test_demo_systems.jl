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

end
