@testset "test_transferfunction" begin
# Naming convention:
# ------------------
# {type}_{dims}
# type: C: Continuous, D: Discrete
# dims: "npnuny" (np = # poles)

# CONTINUOUS
C_111 = tf([1, 2], [1, 5])
C_211 = tf([1, 2, 3], [1, 8, 15])
C_212 = tf(vecarray(2, 1,[1, 2, 3], [1, 2]), vecarray(2, 1, [1, 8, 15], [1, 8, 15]))
C_221 = tf(vecarray(1, 2,[1, 2, 3], [1, 2]), vecarray(1, 2, [1, 8, 15], [1, 8, 15]))
C_222 = [C_221; C_221]
C_022 = tf(4*[1 0; 0 1])
s = tf("s")

# DISCRETE
D_111 = tf([1, 2], [1, -0.5], 0.005)
D_211 = tf([1, 2, 3], [1, -0.2, -0.15], 0.005)
D_221 = tf(vecarray(1, 2, [1, 2, 3], [1, 2]), vecarray(1, 2, [1, -0.2, -0.15], [1, -0.2, -0.15]), 0.005)
D_222 = [D_221; D_221]
D_022 = tf(4*[1 0; 0 1], 0.005)
z = tf("z", 0.005)

# TESTS
# Constructors
@test s == tf([1, 0], [1])
@test z == tf([1, 0], [1], 0.005)
@test C_022 == tf(vecarray(2, 2, [4], [0], [0], [4]), vecarray(2, 2, [1], [1], [1], [1]))
@test D_022 == tf(vecarray(2, 2, [4], [0], [0], [4]), vecarray(2, 2, [1], [1], [1], [1]), 0.005)
@test C_022 == [tf(4) 0;0 4]
@test C_022 == tf([4 0;0 4])
@test D_022 == tf([4 0;0 4], 0.005)

# Test equality
@test tf([1,2], [2,3,4]) == tf(2*[1,2], 2*[2,3,4])
@test tf([1.0], [2.0,3.0]) == tf(π*[1.0], π*[2.0,3.0])
@test tf([1.0+2.0im], [2.0+im,3.0]) == tf((π+im)*[1+2.0im], (π+im)*[2.0+im,3.0])

# Test approximate equlity
# rtol should just be on the order of ϵ, no particular reason that exactly ϵ
# would work, but apparently it does
ϵ = 1e-14
@test tf([1,2], [2,3,4]) != tf(2*[1,2], (2+ϵ)*[2,3,4])
@test tf([1,2], [2,3,4]) ≈   tf(2*[1,2], (2+ϵ)*[2,3,4]) rtol=ϵ
@test tf([1.0], [2.0,3.0]) != tf((π+ϵ)*[1.0], π*[2.0,3.0])
@test tf([1.0], [2.0,3.0]) ≈ tf((π+ϵ)*[1.0], π*[2.0,3.0]) rtol=ϵ
@test tf([1.0+2.0im], [2.0+im,3.0]) != tf((π+(1+ϵ)im)*[1+2.0im], (π+ϵ+im)*[2.0+(1+ϵ*im),3.0])
@test tf([1.0+2.0im], [2.0+im,3.0]) ≈ tf((π+(1+ϵ)im)*[1+2.0im], (π+ϵ+im)*[2.0+(1+ϵ)*im,3.0]) rtol=ϵ

# Addition
@test C_111 + C_111 == tf([2,14,20], [1,10,25])
@test C_222 + C_222 == tf(vecarray(2, 2, [2,20,68,108,90], [0,2,20,62,60],
              [2,20,68,108,90], [0,2,20,62,60]),
  vecarray(2, 2, [1,16,94,240,225], [1,16,94,240,225],
              [1,16,94,240,225], [1,16,94,240,225]))
@test C_222 + 1 == tf(vecarray(2, 2, [2,10,18], [1,9,17], [2,10,18], [1,9,17]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
@test D_111 + D_111 == tf([2,3,-2], [1,-1,0.25], 0.005)

# Subtraction
@test C_111 - C_211 == tf([0,3,18,15], [1,13,55,75])
@test 1 - C_222 == tf(vecarray(2, 2, [0,6,12], [1,7,13], [0,6,12], [1,7,13]), vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
# We are not doing enough to identify zero numerator here
@test_broken D_111 - D_211 - tf([0,0.3,-2.55,1.2], [1,-0.7,-0.05,0.075], 0.005) == tf([0.0], [1], 0.005)

# Multiplication
@test C_111 * C_221 == tf(vecarray(1, 2, [1,4,7,6], [0,1,4,4]),
  vecarray(1, 2, [1,13,55,75], [1,13,55,75]))
@test C_212 * C_111 == tf(vecarray(2, 1, [1,4,7,6], [0,1,4,4]),
  vecarray(2, 1, [1,13,55,75], [1,13,55,75]))
@test 4*C_222 == tf(vecarray(2, 2, [4,8,12], [0,4,8], [4,8,12], [0,4,8]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
 # We are not doing enough to identify zero numerator here
@test_broken D_111 * D_221 - tf(vecarray(1, 2, [1,4,7,6], [0,1,4,4]),
  vecarray(1, 2, [1,-0.7,-0.05,0.075], [1.0,-0.7,-0.05,0.075]), 0.005) ==
tf(vecarray(1, 2, [0], [0]), vecarray(1, 2, [1], [1]), 0.005)

# Division
@test 1/C_111 == tf([1,5], [1,2])
@test C_212/C_111 == tf(vecarray(2, 1, [1,7,13,15], [0,1,7,10]),
  vecarray(2, 1, [1,10,31,30], [1,10,31,30]))
@test 1/D_111 == tf([1.0,-0.5], [1.0,2.0], 0.005)

# Indexing
@test size(C_222) == (2, 2)
@test size(C_212) == (2, 1)
@test C_222[1,1] == tf([1, 2, 3], [1, 8, 15])
@test C_222[1:1,1] == tf([1, 2, 3], [1, 8, 15])
@test C_222[1,1:2] == C_221
@test size(C_222[1,[]]) == (1,0)

# Accessing Ts through .Ts
@test D_111.Ts == 0.005

# propertynames
@test propertynames(C_111) == (:matrix, :timeevol, :nu, :ny)
@test propertynames(D_111) == (:matrix, :timeevol, :nu, :ny, :Ts)


# Errors
@test_throws ErrorException tf(vecarray(1, 1, [1,7,13,15]),
   vecarray(2, 1, [1,10,31,30], [1,10,31,30]))


# Printing
res = ("TransferFunction{Continuous,ControlSystems.SisoRational{Int64}}\nInput 1 to output 1\ns^2 + 2s + 3\n-------------\ns^2 + 8s + 15\n\nInput 1 to output 2\ns^2 + 2s + 3\n-------------\ns^2 + 8s + 15\n\nInput 2 to output 1\n    s + 2\n-------------\ns^2 + 8s + 15\n\nInput 2 to output 2\n    s + 2\n-------------\ns^2 + 8s + 15\n\nContinuous-time transfer function model")
@test dropwhitespace(sprint(show, C_222)) == dropwhitespace(res)
res = ("TransferFunction{Discrete{Float64},ControlSystems.SisoRational{Float64}}\nInput 1 to output 1\n1.0z^2 + 2.0z + 3.0\n--------------------\n1.0z^2 - 0.2z - 0.15\n\nInput 1 to output 2\n1.0z^2 + 2.0z + 3.0\n--------------------\n1.0z^2 - 0.2z - 0.15\n\nInput 2 to output 1\n     1.0z + 2.0\n--------------------\n1.0z^2 - 0.2z - 0.15\n\nInput 2 to output 2\n     1.0z + 2.0\n--------------------\n1.0z^2 - 0.2z - 0.15\n\nSample Time: 0.005 (seconds)\nDiscrete-time transfer function model")
@test dropwhitespace(sprint(show, D_222)) == dropwhitespace(res)


@test tf(zpk([1.0 2; 3 4])) == tf([1 2; 3 4])


# Type stability Continuous-time
@test eltype(fill(tf("s"),2)) <: TransferFunction
@test eltype(fill(tf([1],[1,1]),2)) <: TransferFunction
@test eltype(fill(tf(1.0,[1,1]),2)) <: TransferFunction
@test eltype(fill(tf([1 2; 3 4]),2)) <: TransferFunction
@test eltype(fill(tf(1)+tf(2),2)) <: TransferFunction
@test eltype(fill(tf(1)/tf(2),2)) <: TransferFunction
@test eltype(fill(tf(1)+1,2)) <: TransferFunction

# Type stability Discrete-time
@test eltype(fill(tf("z",1.0),2)) <: TransferFunction
@test eltype(fill(tf([1],[1,1],1),2)) <: TransferFunction
@test eltype(fill(tf(1.0,[1,1],1),2)) <: TransferFunction
@test eltype(fill(tf([1 2; 3 4],1),2)) <: TransferFunction
@test eltype(fill(tf(1,1)+tf(2,1),2)) <: TransferFunction
@test eltype(fill(tf(1,1)/tf(2,1),2)) <: TransferFunction
@test eltype(fill(tf(1,1.0)+1,2)) <: TransferFunction

# Errors
@test_throws ErrorException C_111 + C_222             # Dimension mismatch
@test_throws ErrorException C_111 - C_222             # Dimension mismatch
@test_throws ErrorException C_111 * C_222             # Dimension mismatch
@test_throws ErrorException [s 0; 1]                  # Dimension mismatch
@test_throws ErrorException D_111 + C_111             # Sampling time mismatch
@test_throws ErrorException D_111 - C_111             # Sampling time mismatch
@test_throws ErrorException D_111 * C_111             # Sampling time mismatch
D_diffTs = tf([1], [2], 0.1)
@test_throws ErrorException D_111 + D_diffTs          # Sampling time mismatch
@test_throws ErrorException D_111 - D_diffTs          # Sampling time mismatch
@test_throws ErrorException D_111 * D_diffTs          # Sampling time mismatch
@test_throws ErrorException tf([1], [2], -0.1)        # Negative samping time
@test_throws ErrorException tf("s", 0.01)             # s creation can't be discrete
@test_throws ErrorException tf("z", 0)                # z creation can't be continuous
@test_throws ErrorException tf("z")                   # z creation can't be continuous

@test [z 0] == [tf("z", 0.005) tf(0, 0.005)]

# Test polynomial inputs
Polynomial = ControlSystems.Polynomials.Polynomial
v1, v2, v3, v4 = [1,2,3], [4,5,6], [7,8], [9,10,11]
p1, p2, p3, p4 = Polynomial.(reverse.((v1,v2,v3,v4)))

S1v = tf(v1, v2)
S2v = tf(v3, v4)
Sv = [S1v S2v; 2S1v 3S2v]

# Test polynomial SISO constructor
S1p = tf(p1,p2)
# Test polynomial matrix MIMO constructor
Sp = tf([p1 p3; 2p1 3p3], [p2 p4; p2 p4])

## Test specifying time with TimeEvolution struct
a = [1.0, 2, 1.0]
b = [1.0, 3.0]
@test tf(b, a) == tf(Polynomial(reverse(b)), Polynomial(reverse(a)))
@test tf(b, a) == tf(Polynomial(reverse(b)), Polynomial(reverse(a)), Continuous())
@test tf(b, a, 1.5) == tf(Polynomial(reverse(b)), Polynomial(reverse(a)), 1.5)
@test tf(b, a, 1.5) == tf(Polynomial(reverse(b)), Polynomial(reverse(a)), Discrete(1.5))




@test S1v == S1p
@test Sv == Sp

# Test constructors
@test ControlSystems.SisoRational{Float64}(v1,v2) isa ControlSystems.SisoRational{Float64}
@test ControlSystems.SisoRational{Float64}(v1,1.0*v2) isa ControlSystems.SisoRational{Float64}
@test ControlSystems.SisoRational{Int64}(v1,1.0*v2) isa ControlSystems.SisoRational{Int64}

@test_throws InexactError ControlSystems.SisoRational{Int64}(v1,1.5*v2)

end
