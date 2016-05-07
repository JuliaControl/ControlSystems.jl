module TestTransferFunction
using CustomTest
using ControlSystems

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
C_022 = tf(4*eye(2))
s = tf("s")

# DISCRETE
D_111 = tf([1, 2], [1, -0.5], 0.005)
D_211 = tf([1, 2, 3], [1, -0.2, -0.15], 0.005)
D_221 = tf(vecarray(1, 2, [1, 2, 3], [1, 2]), vecarray(1, 2, [1, -0.2, -0.15], [1, -0.2, -0.15]), 0.005)
D_222 = [D_221; D_221]
D_022 = tf(4*eye(2), 0.005)
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

# Addition
@test C_111 + C_111 == tf([2,14,20], [1,10,25])
@test C_222 + C_222 ==
tf(vecarray(2, 2, [2,20,68,108,90], [0,2,20,62,60],
              [2,20,68,108,90], [0,2,20,62,60]),
  vecarray(2, 2, [1,16,94,240,225], [1,16,94,240,225],
              [1,16,94,240,225], [1,16,94,240,225]))
@test C_222 + 1 == tf(vecarray(2, 2, [2,10,18], [1,9,17], [2,10,18], [1,9,17]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
@test D_111 + D_111 == tf([2,3,-2], [1,-1,0.25], 0.005)

# Subtraction
@test C_111 - C_211 == tf([0,3,18,15], [1,13,55,75])
@test 1 - C_222 == tf(vecarray(2, 2, [0,6,12], [1,7,13], [0,6,12], [1,7,13]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
@test D_111 - D_211 - tf([0,0.3,-2.55,1.2], [1,-0.7,-0.05,0.075], 0.005) ==
    tf([0], [1], 0.005)

# Multiplication
@test C_111 * C_221 == tf(vecarray(1, 2, [1,4,7,6], [0,1,4,4]),
  vecarray(1, 2, [1,13,55,75], [1,13,55,75]))
@test C_212 * C_111 == tf(vecarray(2, 1, [1,4,7,6], [0,1,4,4]),
  vecarray(2, 1, [1,13,55,75], [1,13,55,75]))
@test 4*C_222 == tf(vecarray(2, 2, [4,8,12], [0,4,8], [4,8,12], [0,4,8]),
  vecarray(2, 2, [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
@test D_111 * D_221 - tf(vecarray(1, 2, [1,4,7,6], [0,1,4,4]),
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

# Printing
res = ("TransferFunction:\nInput 1 to Output 1\ns^2 + 2.0s + 3.0\n-----------"*
       "------\ns^2 + 8.0s + 15.0\n\nInput 1 to Output 2\ns^2 + 2.0s + 3.0\n-"*
       "----------------\ns^2 + 8.0s + 15.0\n\nInput 2 to Output 1\n     s + "*
       "2.0\n-----------------\ns^2 + 8.0s + 15.0\n\nInput 2 to Output 2\n   "*
       "  s + 2.0\n-----------------\ns^2 + 8.0s + 15.0\n\nContinuous-time"*
       " transfer function model")
@test sprint(show, C_222) == res
res = ("TransferFunction:\nInput 1 to Output 1\nz^2 + 2.0z + 3.0\n-----------"*
       "------\nz^2 - 0.2z - 0.15\n\nInput 1 to Output 2\nz^2 + 2.0z + 3.0\n-"*
       "----------------\nz^2 - 0.2z - 0.15\n\nInput 2 to Output 1\n     z + "*
       "2.0\n-----------------\nz^2 - 0.2z - 0.15\n\nInput 2 to Output 2\n   "*
       "  z + 2.0\n-----------------\nz^2 - 0.2z - 0.15\n\nSample Time: 0.005 "*
       "(seconds)\nDiscrete-time transfer function model")
@test sprint(show, D_222) == res

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
@test_err C_111 + C_222             # Dimension mismatch
@test_err C_111 - C_222             # Dimension mismatch
@test_err C_111 * C_222             # Dimension mismatch
@test_err [s 0; 1]                  # Dimension mismatch
@test_err D_111 + C_111             # Sampling time mismatch
@test_err D_111 - C_111             # Sampling time mismatch
@test_err D_111 * C_111             # Sampling time mismatch
D_diffTs = tf([1], [2], 0.1)
@test_err D_111 + D_diffTs          # Sampling time mismatch
@test_err D_111 - D_diffTs          # Sampling time mismatch
@test_err D_111 * D_diffTs          # Sampling time mismatch
@test_err tf([1], [2], -0.1)        # Negative samping time
@test_err tf("s", 0.01)             # s creation can't be discrete
@test_err tf("z", 0)                # z creation can't be continuous
@test_err tf("z")                   # z creation can't be continuous
@test_err [z 0]                     # Sampling time mismatch (inferec could be implemented)

# Test numpoly, numvec, denpoly, denvec for SisoRational

# Test deprecation (is not an error)
#@test_err num(C_111.matrix[1,1])
#@test_err den(C_111.matrix[1,1])

@test numvec(C_111.matrix[1,1]) == [1, 2]
@test denvec(C_111.matrix[1,1]) == [1, 5]

vecs = Array{Array{Float64,1},2}(1,2)
vecs[1,1] = [1,2,3]; vecs[1,2] = [1, 2]
@test numvec(D_221) == vecs
vecs[1,1] = [1,-0.2,-0.15]; vecs[1,2] = [1, -0.2, -0.15]
@test denvec(D_221) == vecs

@test size(numpoly(D_221)) == (1,2)
@test numpoly(D_221) == [ControlSystems.Poly([1, 2, 3.0]) ControlSystems.Poly([1, 2.0])]
@test size(denpoly(D_221)) == (1,2)
@test denpoly(D_221) == [ControlSystems.Poly([1, -0.2, -0.15]) ControlSystems.Poly([1, -0.2, -0.15])]

#After switch to polynomials
#@test numpoly(D_221) == [Polynomials.Poly([3.0, 2, 1]) Polynomials.Poly([2.0, 1])])
#@test denpoly(D_221) == [Polynomials.Poly([-0.15, -0.2, 1]) Polynomials.Poly([-0.15, -0.2, 1])]
end
