using ControlSystems
using ControlSystems: to_static, to_sized, freqresp_nohess, freqresp_nohess!
using StaticArrays



@test to_static(1) == 1
@test to_static([1,2]) === SMatrix{2,1,Int,2}(1,2)
@test to_static(SMatrix{2,1,Int,2}(1,2)) === SMatrix{2,1,Int,2}(1,2)
@test to_sized(1) == 1



sys = ss([-5 0; 0 -5],[2; 2],[3 3],[0])
ssys = StaticStateSpace(sys)

@test_nowarn display(ssys)

@test ssys.A isa SMatrix{2,2,eltype(sys.A),4}
@test ssys.B isa SMatrix{2,1,eltype(sys.A),2}
@test ssys.C isa SMatrix{1,2,eltype(sys.A),2}
@test ssys.D isa SMatrix{1,1,eltype(sys.A),1}

@test ssys == sys


ssys = to_sized(sys)
@test_nowarn display(ssys)


@test ssys.A isa SizedMatrix{2,2,eltype(sys.A),2}
@test ssys.B isa SizedMatrix{2,1,eltype(sys.A),2}
@test ssys.C isa SizedMatrix{1,2,eltype(sys.A),2}
@test ssys.D isa SizedMatrix{1,1,eltype(sys.A),2}

@test ssys == sys


sys = ss([-5 0; 0 -5],[2; 2],[3 3],[0], 0.1)
ssys = StaticStateSpace(sys)
@test ssys.Ts == 0.1


@test ssys.A isa SMatrix{2,2,eltype(sys.A),4}
@test ssys.B isa SMatrix{2,1,eltype(sys.A),2}
@test ssys.C isa SMatrix{1,2,eltype(sys.A),2}
@test ssys.D isa SMatrix{1,1,eltype(sys.A),1}

@test ssys == sys


sys = ss(1)
ssys = StaticStateSpace(sys)

@test ssys.A isa SMatrix{0,0,eltype(sys.A),0}
@test ssys.B isa SMatrix{0,1,eltype(sys.A),0}
@test ssys.C isa SMatrix{1,0,eltype(sys.A),0}
@test ssys.D isa SMatrix{1,1,eltype(sys.A),1}

@test ssys == sys


## TransferFunction
G = tf(1, [1,1])
sG = StaticStateSpace(G)
@test sG isa HeteroStateSpace{Continuous, <:SMatrix}
@test sG == ss(G)

@test G * sG == ss(G*ss(G)) == ss(G) * sG
promote_type(typeof(ss(G)), typeof(sG)) == StateSpace{Continuous, Float64}
promote_type(typeof(G), typeof(sG)) == StateSpace{Continuous, Float64}

## Delay system

# sys = DemoSystems.woodberry()
# ssys = StaticStateSpace(sys)


## Multiplication

sys1 = StaticStateSpace(ssrand(1,1,5,proper=false))
sys2 = StaticStateSpace(ssrand(1,1,5,proper=true))

s1 = ssrand(1,1,5,proper=false)
s2 = HeteroStateSpace(ssrand(1,1,5,proper=true))


sysbig = ssrand(Float32, 1, 1, 1)



static_systems = Any[sys1, sys2]
systems = Any[sys1, sys2, s1, s2, sysbig]

for sys1 in systems, sys2 in systems
    # @show typeof.((sys1, sys2))
    @test sys1*sys2 == ss(sys1)*ss(sys2)
    @test sys2*sys1 == ss(sys2)*ss(sys1)
end


for sys1 in static_systems, sys2 in static_systems
    # @show typeof.((sys1, sys2))
    @inferred sys1*sys2
    @inferred sys2*sys1
end


sys12 = sys1*sys2
@test sys12.A isa SMatrix
@test sys12.B isa SMatrix
@test sys12.C isa SMatrix
@test sys12.D isa SMatrix



## Feedback
sys = StaticStateSpace(ssrand(1,1,2,proper=true))

eye = SMatrix{1,1,Float64,1}(1.0)

@inferred feedback(HeteroStateSpace(eye), sys)
@inferred feedback(sys, HeteroStateSpace(eye))
@inferred feedback(sys, sys)

@test feedback(HeteroStateSpace(eye), sys) == feedback(1, ss(sys))

eye, sys2 = promote(1, sys)
@test eye === HeteroStateSpace(eye)
@test sys2 === sys
@inferred promote(sys, 1)
@inferred feedback(1, sys)

@test promote_type(typeof(sys), typeof(1)) == HeteroStateSpace{Continuous}
promote_type(typeof(StaticStateSpace(ssrand(1,1,1,Ts=0.1))), typeof(1)) == HeteroStateSpace{Discrete{Float64}}

sys = StaticStateSpace(ssrand(1,1,2,proper=false))
@test feedback(sys,sys) ≈ feedback(ss(sys),ss(sys))
@inferred  feedback(sys,sys)

## Freqresp

w = exp10.(LinRange(-2, 2, 20000))
R = zeros(ComplexF64, 1,1, length(w))
syss = StaticStateSpace(ssrand(1,1,5,proper=false))
sys = ss(syss)

freqresp_nohess!(R, syss, w) # precompile
@test @allocated(freqresp_nohess!(R, syss, w)) == 0

@test freqresp_nohess!(R, syss, w) ≈ freqresp(ss(syss), w)

D = StaticStateSpace(randn(2,2))
@test freqresp_nohess(D, w) ≈ freqresp_nohess(ss(D), w)

# Benchmarks with length(w) = 20000
# @btime freqresp_nohess!(R, syss, w);
# 29.660 ms (240001 allocations: 409.55 MiB)
# 27.232 ms (200001 allocations: 408.02 MiB) # loop A
# 7.018 ms (0 allocations: 0 bytes) # Separate Static method
# 6.917 ms (0 allocations: 0 bytes) Float64 isntead of Complex C
# 638.585 μs (73 allocations: 9.70 KiB) # @threads
# 621.413 μs (0 allocations: 0 bytes) # Polyester.@batch

# @btime freqresp_nohess!(R, sys, w);
# 15.132 ms (240002 allocations: 22.89 MiB)
# 13.968 ms (160002 allocations: 19.84 MiB) # loop over Adiag