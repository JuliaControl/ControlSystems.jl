@testset "test_zpk" begin
# Naming convention:
# ------------------
# {type}_{dims}
# type: C: Continuous, D: Discrete
# dims: "npnuny" (np = # poles)

# CONTINUOUS
C_111 = zpk([-2],[-5],1)
C_211 = zpk([-1+sqrt(2)im,-1-sqrt(2)im], [-5,-3],1)
C_212 = zpk(vecarray(2, 1, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(2, 1, [-5,-3], [-5,-3]), fill(1, 2, 1))
C_221 = zpk(vecarray(1, 2, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(1, 2, [-5,-3], [-5,-3]), fill(1, 1, 2))
C_222 = [C_221; C_221]
C_022 = zpk([4 0; 0 4])
s = zpk("s")

# DISCRETE
D_111 = zpk([-2], [0.5], 1, 0.005)
D_211 = zpk([-1+sqrt(2)im,-1-sqrt(2)im], [-0.3, 0.5], 1, 0.005)
D_221 = zpk(vecarray(1, 2, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(1, 2, [-0.3, 0.5], [-0.3, 0.5]), [1 1], 0.005)
D_222 = [D_221; D_221]
D_022 = zpk([4 0; 0 4], 0.005)
z = zpk("z", 0.005)

# TESTS
# Constructors
@test s == zpk([0], Int64[], 1)
@test z == zpk([0], Int64[], 1, 0.005)
@test C_022 == zpk(vecarray(Int64, 2, 2, Int64[], Int64[], Int64[], Int64[]), vecarray(Int64, 2, 2, Int64[], Int64[], Int64[], Int64[]), [4 0; 0 4])
@test D_022 == zpk(vecarray(2, 2, Int64[], Int64[], Int64[], Int64[]), vecarray(2, 2, Int64[], Int64[], Int64[], Int64[]), [4 0; 0 4], 0.005)
@test C_022 == [zpk(4) 0;0 4]

@test D_022 == [zpk(4, 0.005) 0; 0 4]

# Test constructors with empty matrices of type Any
@test zpk([-1.0], [], 1.0) == zpk([-1.0], Float64[], 1.0)
@test zpk([], [-1.0], 1.0) == zpk(Float64[], [-1.0], 1.0)
@test zpk([], [], 1.0) == zpk(Float64[], Float64[], 1.0)
@test zpk([], [1.0+im,1.0-im], 1.0) == zpk(ComplexF64[], [1.0+im,1.0-im], 1.0)

# Test constructors with different conjugate ordering
@test zpk([], [1.0-im,1.0+im], 1.0) == zpk(ComplexF64[], [1.0+im,1.0-im], 1.0)
@test zpk([], [1.0+im,1.0,1.0-im], 1.0) == zpk(ComplexF64[], [1.0+im,1.0-im,1.0], 1.0)
@test zpk([1.0-im,1.0,1.0+im], [], 1.0) == zpk([1.0+im,1.0-im,1.0], ComplexF64[], 1.0)
@test zpk([], [1.0+im,1.0,1.0+im,1.0-im,1.0-im], 1.0) == zpk(ComplexF64[], [1.0+im,1.0-im,1.0+im,1.0-im,1.0], 1.0)
@test_throws AssertionError zpk([], [1.01+im,1.0-im], 1.0) 


#TODO improve polynomial accuracy se these are equal
# Addition
@test C_111 + C_111 ≈ zpk([-2,-5], [-5,-5], 2)
@test C_222 + C_222 ≈
  zpk(vecarray(2, 2, [-1+sqrt(2)im,-1-sqrt(2)im], [-2], [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(2, 2, [-5,-3], [-5,-3], [-5,-3], [-5,-3]), [2 2;2 2])
z1 = [-2.5+sqrt(9-2.5^2)im, -2.5-sqrt(9-2.5^2)im]
z2 = [-4.5-sqrt(4.5^2-17), -4.5+sqrt(4.5^2-17)]


@test C_222 + 1 ≈ zpk(vecarray(2, 2, z1, z2, z1, z2), vecarray(2, 2, [-3, -5], [-3, -5], [-3, -5], [-3, -5]), [2 1; 2 1])
@test D_111 + D_111 ≈ zpk([-2],[0.5],2,0.005)

# Subtraction
@test C_111 - C_211 ≈ zpk([-1], [-5,-3.0], 3.0)
@test minreal(zpk(tf([0,3.0,18,15], [1,13,55,75])), 1e-6) ≈ zpk([-1.0], [-5,-3.0], 3.0)

@test 1 - C_222 ≈ zpk(tf(vecarray(2, 2, [0,6,12], [1,7,13], [0,6,12], [1,7,13]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15])))
@test D_111 - D_211 - zpk(tf([0,0.3,-2.55,1.2], [1,-0.7,-0.05,0.075], 0.005)) ≈
  zpk(tf([0], [1], 0.005))

# Multiplication
@test C_111 * C_221 ≈ zpk(tf(vecarray(1, 2, [1,4,7,6], [0,1,4,4]),
  vecarray(1, 2, [1,13,55,75], [1,13,55,75])))
@test C_212 * C_111 ≈ zpk(tf(vecarray(2, 1, [1,4,7,6], [0,1,4,4]),
  vecarray(2, 1, [1,13,55,75], [1,13,55,75])))
@test 4*C_222 ≈ zpk(tf(vecarray(2, 2, [4,8,12], [0,4,8], [4,8,12], [0,4,8]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15])))
@test D_111 * D_221 - zpk(tf(vecarray(1, 2, [1,4,7,6], [0,1,4,4]),
  vecarray(1, 2, [1,-0.7,-0.05,0.075], [1.0,-0.7,-0.05,0.075]), 0.005)) ≈
  zpk(tf(vecarray(1, 2, [0], [0]), vecarray(1, 2, [1], [1]), 0.005))

# Division
@test 1/C_111 ≈ zpk(tf([1,5], [1,2]))
@test C_212/C_111 ≈ zpk(tf(vecarray(2, 1, [1,7,13,15], [0,1,7,10]),
  vecarray(2, 1, [1,10,31,30], [1,10,31,30])))
@test 1/D_111 ≈ zpk(tf([1.0,-0.5], [1.0,2.0], 0.005))

# Indexing
@test size(C_222) == (2, 2)
@test size(C_212) == (2, 1)
@test C_222[1,1] ≈ zpk(tf([1, 2, 3], [1, 8, 15]))
@test C_222[1:1,1] ≈ zpk(tf([1, 2, 3], [1, 8, 15]))
@test C_222[1,1:2] == C_221
@test size(C_222[1,Int64[]]) == (1,0)

# Test that number of poles matter
@test !(zpk(Int64[],[1,1],1) == zpk(Int64[],[1],1))

# Test SispZpk and SisoTf operations
C_212_tf = tf(C_212)
C_111_tf = tf(C_111)
# Add
@test C_212_tf + C_212 ≈ 2*C_212
@test C_212 + C_212_tf ≈ 2*C_212
# Approx
@test C_212_tf + C_212 ≈ 2*C_212_tf
# Minus
@test 2*C_212_tf - C_212 ≈ C_212
# Multiply
@test C_212_tf*C_111 ≈ C_212*C_111
@test C_212_tf*C_111 ≈ C_212_tf*C_111_tf



## Test specifying time with TimeEvolution struct
zvec = [0.5]
pvec = [-0.2 + 0.3im, -0.2 - 0.3im]
k = 0.3
@test zpk(zvec, pvec, k) == zpk(zvec, pvec, k, Continuous())
@test zpk(zvec, pvec, k, 0.2) == zpk(zvec, pvec, k, Discrete(0.2))



# TODO test printing when it is implemented better

# Tests of minreal
@test minreal(zpk([-5.0], [-5.0, -5.0], 1.0)) == zpk(Float64[], [-5.0], 1.0)
@test minreal(zpk([-1.0, -2.0], [-2.0, -2.0], 1.0)) == zpk([-1.0], [-2.0], 1.0)

# Tests of minreal (roudning and cancelling complex roots with real roots)
@test minreal(zpk([-2.0+1e-10im,-2.0-1e-10im], [-2.0], 1.0)) == zpk([-2.0], Float64[], 1.0)
@test minreal(zpk([-2.0], [-2.0+1e-10im,-2.0-1e-10im], 1.0)) == zpk(Float64[], [-2.0], 1.0)
@test minreal(zpk([-1.0, -2.0, -2.0+1e-10im,-2.0-1e-10im], [-2.0, -2.0], 1.0)) == zpk([-1.0, -2.0], Float64[], 1.0)
@test minreal(zpk([-1.0, -2.0+1e-10im,-2.0-1e-10im,-2.0], [-2.0, -2.0], 1.0)) == zpk([-1.0, -2.0], Float64[], 1.0)
@test minreal(zpk([-2.0, -2.0], [-1.0, -2.0, -2.0+1e-10im,-2.0-1e-10im], 1.0)) ==  zpk(Float64[], [-1.0, -2.0], 1.0)
@test minreal(zpk([-2.0, -2.0], [-1.0, -2.0+1e-10im,-2.0-1e-10im,-2.0], 1.0)) == zpk(Float64[], [-1.0, -2.0], 1.0)

# Test of minreal (systems with no poles / no zeros)
@test minreal(zpk([-1.0, -2.0], Float64[], 2.5)) == zpk([-1.0, -2.0], Float64[], 2.5)
@test minreal(zpk(Float64[], [-1.0, -2.0], 2.5)) == zpk(Float64[], [-1.0, -2.0], 2.5)

# Test type inference
@test eltype(fill(zpk("s"),2)) <: TransferFunction
@test eltype(fill(zpk([1],[1,1],1),2)) <: TransferFunction
@test eltype(fill(zpk(Int64[],[1,1],1.0),2)) <: TransferFunction
@test eltype(fill(zpk([1 2; 3 4],0.001),2)) <: TransferFunction
@test eltype(fill(zpk(1)+zpk(2),2)) <: TransferFunction
@test eltype(fill(zpk(1,0.005)/zpk(2, 0.005),2)) <: TransferFunction
@test eltype(fill(zpk(1)+1,2)) <: TransferFunction

@test eltype([tf(1,1), zpk(1,1)]) <: TransferFunction

zpk(tf([1 2; 3 4])) == zpk([1 2; 3 4])

# Errors
@test_throws ErrorException C_111 + C_222             # Dimension mismatch
@test_throws ErrorException C_111 - C_222             # Dimension mismatch
@test_throws ErrorException C_111 * C_222             # Dimension mismatch
@test_throws ErrorException [s 0; 1]                  # Dimension mismatch
@test_throws ErrorException D_111 + C_111             # Sampling time mismatch
@test_throws ErrorException D_111 - C_111             # Sampling time mismatch
@test_throws ErrorException D_111 * C_111             # Sampling time mismatch
D_diffTs = zpk(tf([1], [2], 0.1))
@test_throws ErrorException D_111 + D_diffTs          # Sampling time mismatch
@test_throws ErrorException D_111 - D_diffTs          # Sampling time mismatch
@test_throws ErrorException D_111 * D_diffTs          # Sampling time mismatch
@test_throws ErrorException zpk([1], [2], 1, -0.1)        # Negative samping time
@test_throws ErrorException zpk("s", 0.01)             # s creation can't be discrete
@test_throws ErrorException zpk("z", 0)                # z creation can't be continuous
@test_throws ErrorException zpk("z")                   # z creation can't be continuous

@test [z 0] == [zpk([0], Int64[], 1, 0.005) zpk([], [], 0, 0.005)]


@test typeof(zpk(tf([1], [2], 0.1))) == TransferFunction{Discrete{Float64},ControlSystems.SisoZpk{Float64,Complex{Float64}}}
@test typeof(zpk([-0.5], [], 1)) == TransferFunction{Continuous,ControlSystems.SisoZpk{Float64,Float64}}
end
