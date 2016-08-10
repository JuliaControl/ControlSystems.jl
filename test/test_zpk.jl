module TestZpk
using CustomTest
using Base.Test
using ControlSystems

# Naming convention:
# ------------------
# {type}_{dims}
# type: C: Continuous, D: Discrete
# dims: "npnuny" (np = # poles)

# CONTINUOUS
C_111 = zpk([-2],[-5],1)
C_211 = zpk([-1+sqrt(2)im,-1-sqrt(2)im], [-5,-3],1)
C_212 = zpk(vecarray(2, 1, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(2, 1, [-5,-3], [-5,-3]), [1;1])
C_221 = zpk(vecarray(1, 2, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(1, 2, [-5,-3], [-5,-3]), [1 1])
C_222 = [C_221; C_221]
C_022 = zpk(4eye(2))
s = zpk("s")

# DISCRETE
D_111 = zpk([-2], [0.5], 1, 0.005)
D_211 = zpk([-1+sqrt(2)im,-1-sqrt(2)im], [-0.3, 0.5], 1, 0.005)
D_221 = zpk(vecarray(1, 2, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(1, 2, [-0.3, 0.5], [-0.3, 0.5]), [1 1], 0.005)
D_222 = [D_221; D_221]
D_022 = zpk(4eye(2), 0.005)
z = zpk("z", 0.005)

# TESTS
# Constructors
@test s == zpk([0], [], 1)
@test z == zpk([0], [], 1, 0.005)
@test C_022 == zpk(vecarray(Int64, 2, 2, [], [], [], []), vecarray(Int64, 2, 2, [], [], [], []), [4 0; 0 4])
@test D_022 == zpk(vecarray(2, 2, [], [], [], []), vecarray(2, 2, [], [], [], []), [4 0; 0 4], 0.005)
@test C_022 == [zpk(4) 0;0 4]
#TODO We might want to fix this
#@test D_022 == [zpk(4, 0.005) 0;0 4])

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
@test C_111 - C_211 ≈ zpk(tf([0,3,18,15], [1,13,55,75]))
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
@test C_222[1,1] == zpk(tf([1, 2, 3], [1, 8, 15]))
@test C_222[1:1,1] == zpk(tf([1, 2, 3], [1, 8, 15]))
@test C_222[1,1:2] == C_221
@test size(C_222[1,[]]) == (1,0)

# Test that number of poles matter
@test !(zpk([],[1,1],1) == zpk([],[1],1))
# TODO test printing when it is implemented better

# Type stability Continuous and discrete time
@test eltype(fill(zpk("s"),2)) <: TransferFunction
@test eltype(fill(zpk([1],[1,1],1),2)) <: TransferFunction
@test eltype(fill(zpk([],[1,1],1.0),2)) <: TransferFunction
@test eltype(fill(zpk([1 2; 3 4],0.001),2)) <: TransferFunction
@test eltype(fill(zpk(1)+zpk(2),2)) <: TransferFunction
@test eltype(fill(zpk(1,0.005)/zpk(2, 0.005),2)) <: TransferFunction
@test eltype(fill(zpk(1)+1,2)) <: TransferFunction

@test eltype([tf(1,1), zpk(1,1)]) <: TransferFunction

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
# Remove this when inferec is implemented
@test_throws ErrorException [z 0]                     # Sampling time mismatch (inferec could be implemented)
end
