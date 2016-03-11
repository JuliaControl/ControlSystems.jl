module TestConnections
using CustomTest
using Control

## STATE SPACE ##
# CONTINUOUS
C_111 = ss([1], [2], [3], [4])
C_211 = ss(eye(2), [1; 2], [1 0], [0])
C_212 = ss(eye(2), [1; 2], eye(2), [0; 0])
C_221 = ss(eye(2), [1 0; 0 2], [1 0], [0 0])
C_222 = ss(eye(2), [1 0; 0 2], eye(2), zeros(2,2))
C_022 = ss(4*eye(2))

# DISCRETE
D_111 = ss([1], [2], [3], [4], 0.005)
D_211 = ss(eye(2), [1; 2], [1 0], [0], 0.005)
D_212 = ss(eye(2), [1;2], eye(2), [0; 0], 0.005)
D_221 = ss(eye(2), [1 0; 0 2], [1 0], [0 0], 0.005)
D_222 = ss(eye(2), [1 0; 0 2], eye(2), zeros(2,2), 0.005)
D_022 = ss(4*eye(2), 0.005)

@test [C_111 C_221] == ss(eye(3), [2 0 0; 0 1 0; 0 0 2], [3 1 0], [4 0 0])
@test [C_111; C_212] == ss(eye(3), [2; 1; 2], [3 0 0; 0 1 0; 0 0 1], [4; 0; 0])
@test append(C_111, C_211) == ss(eye(3), [2 0; 0 1; 0 2], [3 0 0; 0 1 0], [4 0; 0 0])
@test [C_022 C_222] == ss(eye(2), [0 0 1 0; 0 0 0 2], [1 0; 0 1], [4 0 0 0; 0 4 0 0])
@test [C_022; C_222] == ss(eye(2), [1 0; 0 2], [0 0; 0 0; 1 0; 0 1], [4 0; 0 4; 0 0; 0 0])

@test [D_111 D_221] == ss(eye(3), [2 0 0; 0 1 0; 0 0 2], [3 1 0], [4 0 0], 0.005)
@test [D_111; D_212] == ss(eye(3), [2; 1; 2], [3 0 0; 0 1 0; 0 0 1], [4; 0; 0], 0.005)
@test append(D_111, D_211) == ss(eye(3), [2 0; 0 1; 0 2], [3 0 0; 0 1 0], [4 0; 0 0], 0.005)
@test [D_022 D_222] == ss(eye(2), [0 0 1 0; 0 0 0 2], [1 0; 0 1], [4 0 0 0; 0 4 0 0], 0.005)
@test [D_022; D_222] == ss(eye(2), [1 0; 0 2], [0 0; 0 0; 1 0; 0 1], [4 0; 0 4; 0 0; 0 0], 0.005)

@test series(C_111, C_212) == C_212*C_111
@test parallel(C_111, C_211) == C_111 + C_211

# Errors
macro test_err(ex)
    :(@test_throws ErrorException $ex)
end
@test_err [C_111 D_111]                 # Sampling time mismatch
@test_err [C_111; D_111]                # Sampling time mismatch
@test_err append(C_111, D_111)          # Sampling time mismatch
@test_err [C_111 C_212]                 # Dimension mismatch
@test_err [C_111; C_221]                # Dimension mismatch

## TRANSFER FUNCTION ##
# CONTINUOUS
Ctf_111 = tf([1, 2], [1, 5])
Ctf_211 = tf([1, 2, 3], [1, 8, 15])
Ctf_212 = tf(vecarray(2, 1, [1, 2, 3], [1, 2]), vecarray(2, 1, [1, 8, 15], [1, 8, 15]))
Ctf_221 = tf(vecarray(1, 2, [1, 2, 3], [1, 2]), vecarray(1, 2, [1, 8, 15], [1, 8, 15]))
Ctf_222 = [Ctf_221; Ctf_221]
Ctf_022 = tf(4*eye(2))

# DISCRETE
Dtf_111 = tf([1, 2], [1, 5], 0.005)
Dtf_211 = tf([1, 2, 3], [1, 8, 15], 0.005)
Dtf_212 = tf(vecarray(2, 1, [1, 2, 3], [1, 2]), vecarray(2, 1, [1, 8, 15], [1, 8, 15]), 0.005)
Dtf_221 = tf(vecarray(1, 2, [1, 2, 3], [1, 2]), vecarray(1, 2, [1, 8, 15], [1, 8, 15]), 0.005)
Dtf_222 = [Dtf_221; Dtf_221]; Dtf_222.Ts = 0.005
Dtf_022 = tf(4*eye(2), 0.005)

s = tf("s")
@test [Ctf_111 Ctf_221] == tf(vecarray(1, 3, [1,2], [1,2,3], [0,1,2]),
    vecarray(1, 3, [1,5], [1,8,15], [1,8,15]))
@test [Ctf_111; Ctf_212] == tf(vecarray(3, 1, [1,2], [1,2,3], [0,1,2]),
    vecarray(3, 1, [1,5], [1,8,15], [1,8,15]))
@test append(Ctf_111, Ctf_211) == tf(vecarray(2, 2, [1,2], [0], [0], [1,2,3]),
    vecarray(2, 2, [1,5], [1], [1], [1,8,15]));
@test [Ctf_022 Ctf_222] == tf(vecarray(2, 4, [4], [0], [1,2,3], [0,1,2], [0], [4], [1,2,3], [0,1,2]),
    vecarray(2, 4, [1], [1], [1,8,15], [1,8,15], [1], [1], [1,8,15], [1,8,15]))
@test [Ctf_022; Ctf_222] == tf(vecarray(4, 2, [4], [0], [0], [4], [1,2,3], [0,1,2], [1,2,3], [0,1,2]),
    vecarray(4, 2, [1], [1], [1], [1], [1,8,15], [1,8,15], [1,8,15], [1,8,15]))
@test [Ctf_022 Ctf_022] == [[tf(4) 0;0 4] 4*eye(2)]

@test [Dtf_111 Dtf_221] == tf(vecarray(1, 3, [1,2], [1,2,3], [0,1,2]),
    vecarray(1, 3, [1,5], [1,8,15], [1,8,15]), 0.005)
@test [Dtf_111; Dtf_212] == tf(vecarray(3, 1, [1,2], [1,2,3], [0,1,2]),
    vecarray(3, 1, [1,5], [1,8,15], [1,8,15]), 0.005)
@test append(Dtf_111, Dtf_211) == tf(vecarray(2, 2, [1,2], [0], [0], [1,2,3]),
    vecarray(2, 2, [1,5], [1], [1], [1,8,15]), 0.005);
@test [Dtf_022 Dtf_222] == tf(vecarray(2, 4, [4], [0], [1,2,3], [0,1,2], [0], [4], [1,2,3], [0,1,2]),
    vecarray(2, 4, [1], [1], [1,8,15], [1,8,15], [1], [1], [1,8,15], [1,8,15]), 0.005)
@test [Dtf_022; Dtf_222] == tf(vecarray(4, 2, [4], [0], [0], [4], [1,2,3], [0,1,2], [1,2,3], [0,1,2]),
    vecarray(4, 2, [1], [1], [1], [1], [1,8,15], [1,8,15], [1,8,15], [1,8,15]), 0.005)

@test series(Ctf_111, Ctf_212) == tf(vecarray(2, 1, [1,4,7,6], [0,1,4,4]),
    vecarray(2, 1, [1,13,55,75], [1,13,55,75]));
@test parallel(Ctf_111, Ctf_211) == tf([2,17,44,45], [1,13,55,75])

# Combination tf and ss
@test [C_111 Ctf_221] == [C_111 ss(Ctf_221)]
@test [C_111; Ctf_212] == [C_111; ss(Ctf_212)]
@test append(C_111, Ctf_211) == append(C_111, ss(Ctf_211))
@test [D_111 Dtf_221] == [D_111 ss(Dtf_221)]
@test [D_111; Dtf_212] == [D_111; ss(Dtf_212)]
@test append(D_111, Dtf_211) == append(D_111, ss(Dtf_211))
end
