@testset "SimResult" begin
    @info "Testing SimResult"
    
    import ControlSystems: SimResult, ssrand

    y = randn(1,10)
    u = randn(1,10)
    x = randn(1,10)
    t = 1:10
    sys = ssrand(1, 1, 1, Ts=1)

    r = SimResult(y,t,x,u,sys)

    # test getindex
    @test r[1] == y
    @test r[2] == t
    @test r[3] == x
    @test r[4] == u
    @test r[5] == sys

    # test destructuring
    y2,t2,x2,u2 = r
    @test y2 === y
    @test t2 === t
    @test x2 === x
    @test u2 === u

    # test inferability of destructuring
    foo(r::SimResult) = (a,b,c,d) = r
    @inferred foo(r)

    # test properties
    @test r.nx == 1
    @test r.nu == 1
    @test r.ny == 1
end

