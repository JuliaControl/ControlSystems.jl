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

# Printing (default, specific precision, zpk)
res = ("TransferFunction:\n"*
       "Input 1 to Output 1\n"*
       "s^2 + 2*s + 3\n"*
       "--------------\n"*
       "s^2 + 8*s + 15\n"*
       "\n"*
       "Input 1 to Output 2\n"*
       "s^2 + 2*s + 3\n"*
       "--------------\n"*
       "s^2 + 8*s + 15\n"*
       "\n"*
       "Input 2 to Output 1\n"*
       "    s + 2\n"*
       "--------------\n"*
       "s^2 + 8*s + 15\n"*
       "\n"*
       "Input 2 to Output 2\n"*
       "    s + 2\n"*
       "--------------\n"*
       "s^2 + 8*s + 15\n"*
       "\n"*
       "Continuous-time transfer function model")
@test sprint(show, C_222) == res
res = ("TransferFunction:\n"*
       "Input 1 to Output 1\n"*
       "  z^2 + 2*z + 3\n"*
       "------------------\n"*
       "z^2 - 0.2*z - 0.15\n"*
       "\n"*
       "Input 1 to Output 2\n"*
       "  z^2 + 2*z + 3\n"*
       "------------------\n"*
       "z^2 - 0.2*z - 0.15\n"*
       "\n"*
       "Input 2 to Output 1\n"*
       "      z + 2\n"*
       "------------------\n"*
       "z^2 - 0.2*z - 0.15\n"*
       "\n"*
       "Input 2 to Output 2\n"*
       "      z + 2\n"*
       "------------------\n"*
       "z^2 - 0.2*z - 0.15\n"*
       "\n"*
       "Sample Time: 0.005 (seconds)\n"*
       "Discrete-time transfer function model")
@test sprint(show, D_222) == res
res = ("TransferFunction:\n"*
       "Input 1 to Output 1\n"*
       "  s^2 + 2*s + 3\n"*
       "-----------------\n"*
       "s^2 + 8*s + 2e+01\n"*
       "\n"*
       "Input 1 to Output 2\n"*
       "  s^2 + 2*s + 3\n"*
       "-----------------\n"*
       "s^2 + 8*s + 2e+01\n"*
       "\n"*
       "Input 2 to Output 1\n"*
       "      s + 2\n"*
       "-----------------\n"*
       "s^2 + 8*s + 2e+01\n"*
       "\n"*
       "Input 2 to Output 2\n"*
       "      s + 2\n"*
       "-----------------\n"*
       "s^2 + 8*s + 2e+01\n"*
       "\n"*
       "Continuous-time transfer function model")
@test sprint(show, C_222, 1) == res
res = ("TransferFunction:\n"*
       "Input 1 to Output 1\n"*
       "  z^2 + 2*z + 3\n"*
       "-----------------\n"*
       "z^2 - 0.2*z - 0.1\n"*
       "\n"*
       "Input 1 to Output 2\n"*
       "  z^2 + 2*z + 3\n"*
       "-----------------\n"*
       "z^2 - 0.2*z - 0.1\n"*
       "\n"*
       "Input 2 to Output 1\n"*
       "      z + 2\n"*
       "-----------------\n"*
       "z^2 - 0.2*z - 0.1\n"*
       "\n"*
       "Input 2 to Output 2\n"*
       "      z + 2\n"*
       "-----------------\n"*
       "z^2 - 0.2*z - 0.1\n"*
       "\n"*
       "Sample Time: 0.005 (seconds)\n"*
       "Discrete-time transfer function model")
@test sprint(show, D_222, 1) == res
res = ("TransferFunction:\n"*
       "              s - 1.12345\n"*
       "4.45678--------------------------\n"*
       "       (s - 3.67891)(s - 2.54321)\n"*
       "\n"*
       "Continuous-time transfer function model")
@test sprint(show, zpk([1.12345], [2.54321;3.67891], 4.45678)) == res
res = ("TransferFunction:\n"*
       "        s - 1.1\n"*
       "4.5------------------\n"*
       "   (s - 3.7)(s - 2.5)\n"*
       "\n"*
       "Continuous-time transfer function model")
@test sprint(show, zpk([1.12345], [2.54321;3.67891], 4.45678), 2) == res

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
@test_throws ErrorException [z 0]                     # Sampling time mismatch (inferec could be implemented)
end
