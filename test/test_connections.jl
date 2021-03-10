@testset "test_connections" begin
## STATE SPACE ##
# CONTINUOUS
C_111 = ss([1], [2], [3], [4])
C_211 = ss(eye_(2), [1; 2], [1 0], [0])
C_212 = ss(eye_(2), [1; 2], eye_(2), [0; 0])
C_221 = ss(eye_(2), [1 0; 0 2], [1 0], [0 0])
C_222 = ss(eye_(2), [1 0; 0 2], eye_(2), zeros(Int,2,2))
C_022 = ss(4*eye_(2))

# DISCRETE
D_111 = ss([1], [2], [3], [4], 0.005)
D_211 = ss(eye_(2), [1; 2], [1 0], [0], 0.005)
D_212 = ss(eye_(2), [1;2], eye_(2), [0; 0], 0.005)
D_221 = ss(eye_(2), [1 0; 0 2], [1 0], [0 0], 0.005)
D_222 = ss(eye_(2), [1 0; 0 2], eye_(2), zeros(Int,2,2), 0.005)
D_022 = ss(4*eye_(2), 0.005)

@test [C_111 C_221] == ss(eye_(3), [2 0 0; 0 1 0; 0 0 2], [3 1 0], [4 0 0])
@test [C_111; C_212] == ss(eye_(3), [2; 1; 2], [3 0 0; 0 1 0; 0 0 1], [4; 0; 0])
@test append(C_111, C_211) == ss(eye_(3), [2 0; 0 1; 0 2], [3 0 0; 0 1 0], [4 0; 0 0])
@test [C_022 C_222] == ss(eye_(2), [0 0 1 0; 0 0 0 2], [1 0; 0 1], [4 0 0 0; 0 4 0 0])
@test [C_022; C_222] == ss(eye_(2), [1 0; 0 2], [0 0; 0 0; 1 0; 0 1], [4 0; 0 4; 0 0; 0 0])

@test [D_111 D_221] == ss(eye_(3), [2 0 0; 0 1 0; 0 0 2], [3 1 0], [4 0 0], 0.005)
@test [D_111; D_212] == ss(eye_(3), [2; 1; 2], [3 0 0; 0 1 0; 0 0 1], [4; 0; 0], 0.005)
@test append(D_111, D_211) == ss(eye_(3), [2 0; 0 1; 0 2], [3 0 0; 0 1 0], [4 0; 0 0], 0.005)
@test [D_022 D_222] == ss(eye_(2), [0 0 1 0; 0 0 0 2], [1 0; 0 1], [4 0 0 0; 0 4 0 0], 0.005)
@test [D_022; D_222] == ss(eye_(2), [1 0; 0 2], [0 0; 0 0; 1 0; 0 1], [4 0; 0 4; 0 0; 0 0], 0.005)

@test series(C_111, C_212) == C_212*C_111
@test parallel(C_111, C_211) == C_111 + C_211

# Errors
@test_throws ErrorException [C_111 D_111]                 # Sampling time mismatch
@test_throws ErrorException [C_111; D_111]                # Sampling time mismatch
@test_throws ErrorException append(C_111, D_111)          # Sampling time mismatch
@test_throws ErrorException [C_111 C_212]                 # Dimension mismatch
@test_throws ErrorException [C_111; C_221]                # Dimension mismatch

## TRANSFER FUNCTION ##
# CONTINUOUS
Ctf_111 = tf([1, 2], [1, 5])
Ctf_211 = tf([1, 2, 3], [1, 8, 15])
Ctf_212 = tf(vecarray(2, 1, [1, 2, 3], [1, 2]), vecarray(2, 1, [1, 8, 15], [1, 8, 15]))
Ctf_221 = tf(vecarray(1, 2, [1, 2, 3], [1, 2]), vecarray(1, 2, [1, 8, 15], [1, 8, 15]))
Ctf_222 = [Ctf_221; Ctf_221]
Ctf_022 = tf(4*eye_(2))

# DISCRETE
Dtf_111 = tf([1, 2], [1, 5], 0.005)
Dtf_211 = tf([1, 2, 3], [1, 8, 15], 0.005)
Dtf_212 = tf(vecarray(2, 1, [1, 2, 3], [1, 2]), vecarray(2, 1, [1, 8, 15], [1, 8, 15]), 0.005)
Dtf_221 = tf(vecarray(1, 2, [1, 2, 3], [1, 2]), vecarray(1, 2, [1, 8, 15], [1, 8, 15]), 0.005)
Dtf_222 = [Dtf_221; Dtf_221];
Dtf_022 = tf(4*eye_(2), 0.005)

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
@test [Ctf_022 Ctf_022] == [[tf(4) 0;0 4] 4*eye_(2)]

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

# Combination of DelayLtiSystem with TransferFunction and StateSpace
@test [delay(1.0) tf(1, [1,2])] == [delay(1.0) ss(-2.0,1,1,0)]
@test [delay(1.0) zpk([], [-2], 1)] == [delay(1.0) ss(-2.0,1,1,0)]

# hcat and vcat for StateSpace and Matrix
A = [-1.1 -1.2; -1.3 -1.4]
B = [1 2; 3 4]
C = [5 6; 7 8]
D = [1 0; 0 1]
P = ss(A, B, C, D)
@test [P fill(2.5, 2, 1)] == ss(A, [B fill(0, 2, 1)], C, [D fill(2.5, 2, 1)])
@test [fill(2.5, 2, 1) P] == ss(A, [fill(0, 2, 1) B], C, [fill(2.5, 2, 1) D])
@test [P; fill(2.5, 1, 2)] == ss(A, B, [C; fill(0, 1, 2)], [D; fill(2.5, 1, 2)])
@test [fill(2.5, 1, 2); P] == ss(A, B, [fill(0, 1, 2); C], [fill(2.5, 1, 2); D])

# hcat and vcat for StateSpace and Number
P = ss(-1.0, 2.0, 3.0, 4.0)
@test [P 2.5] == ss(-1.0, [2.0 0.0], 3.0, [4.0 2.5])
@test [2.5 P] == ss(-1.0, [0.0 2.0], 3.0, [2.5 4.0])
@test [P; 2.5] == ss(-1.0, 2.0, [3.0; 0.0], [4.0; 2.5])
@test [2.5; P] == ss(-1.0, 2.0, [0.0; 3.0], [2.5; 4.0])

@test [2.5 P 3.5] == ss(-1.0, [0.0 2.0 0.0], 3.0, [2.5 4.0 3.5])
@test [2.5; P; 3.5] == ss(-1.0, 2.0, [0.0; 3.0; 0.0], [2.5; 4.0; 3.5])



# Combination tfRational and sisoZpk
Czpk_111 = zpk([-2],[-5],1)
Czpk_211 = zpk([-1+sqrt(2)im,-1-sqrt(2)im], [-5,-3],1)
Czpk_212 = zpk(vecarray(2, 1, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(2, 1, [-5,-3], [-5,-3]), fill(1, 2, 1))
Czpk_221 = zpk(vecarray(1, 2, [-1+sqrt(2)im,-1-sqrt(2)im], [-2]), vecarray(1, 2, [-5,-3], [-5,-3]), fill(1, 1, 2))
Czpk_222 = [Czpk_221; Czpk_221]
Czpk_022 = [zpk([],[],4) 0; 0 zpk([],[],4)]

#Make sure that we get a vector
arr = Array{typeof(zpk(tf(1))),1}(undef, 2)
arr[1] = zpk(tf(1)); arr[2] = zpk(2);
@test [tf(1), zpk(2)] == arr
arr2 = Array{typeof(tf(1, 0.1)),1}(undef, 2)
arr2[1] = tf(1, 0.1); arr2[2] = tf(2, 0.1);
@test [tf(1, 0.1), tf(2, 0.1)] == arr2
arr3 = Array{typeof(ss(0.)),1}(undef, 3)
arr3[1] = ss(0.); arr3[2] = ss(1.); arr3[3] = ss(2.)
@test [0., zpk(1), ss(2.)] == arr3

arr4 = Array{typeof(ss(0)),1}(undef, 3)
arr4[1] = ss(0); arr4[2] = ss(1); arr4[3] = ss(2)
@test [0., zpk(1), ss(2)] == arr4

@test Czpk_111 ≈ Ctf_111
@test Czpk_211 ≈ Ctf_211
@test Czpk_212 ≈ Ctf_212
@test Czpk_221 ≈ Ctf_221
@test Czpk_222 ≈ Ctf_222
@test Czpk_022 ≈ Ctf_022

@test Czpk_222 ≈ [Ctf_221; Czpk_221]

#This might fail depending on if minreal is used or not
@test (Czpk_211+1) ≈ (Ctf_211+1)


# Concatenation of discrete system with constant
@test [D_111 1.0] == ss([1.0], [2.0 0.0], [3.0], [4.0 1.0], 0.005)
@test [1.0 D_111] == ss([1.0], [0.0 2.0], [3.0], [1.0 4.0], 0.005)
# Type and sample time
@test [D_111 1.0] isa StateSpace{Discrete{Float64},Float64}
@test [D_111 1.0].Ts == 0.005
# Continuous version
@test [C_111 1.0] == ss([1.0], [2.0 0.0], [3.0], [4.0 1.0])
@test [1.0 C_111] == ss([1.0], [0.0 2.0], [3.0], [1.0 4.0])
@test [C_111 1.0] isa StateSpace{Continuous,Float64}
@test [C_111 1.0].Ts == 0.0
@test_logs (:warn,
            "Getting time 0.0 for non-discrete systems is deprecated. Check `isdiscrete` before trying to access time."
            ) [C_111 1.0].Ts
# Concatenation of discrete system with matrix
@test [D_222 fill(1.5, 2, 2)] == [D_222 ss(fill(1.5, 2, 2),0.005)]
@test [C_222 fill(1.5, 2, 2)] == [C_222 ss(fill(1.5, 2, 2))]

# hvcat numbers (second row should be properly handled)
@test [C_111 1.5; 2 3] ==
    [C_111 ss(1.5); ss(2.0) ss(3.0)]
@test [D_111 1.5; 2 3] ==
    [D_111 ss(1.5,0.005); ss(2.0,0.005) ss(3.0,0.005)]
# hvcat matrices
@test [C_222 fill(1.5, 2, 2); fill(2, 2, 2) fill(3, 2, 2)] ==
    [C_222 ss(fill(1.5, 2, 2)); ss(fill(2, 2, 2)) ss(fill(3, 2, 2))]
@test [D_222 fill(1.5, 2, 2); fill(2, 2, 2) fill(3, 2, 2)] ==
    [D_222 ss(fill(1.5, 2, 2),0.005); ss(fill(2, 2, 2),0.005) ss(fill(3, 2, 2),0.005)]


# continuous-time systems
G1 = ss(-9, [2 3], [4; 5], 0)
G2 = ss(-6, 7, 8, 0)
G3 = ss(-1, 1, 1, 1) #  Not strictly proper
K1 = ss(-1, 1, 1, 0)

# discrete-time systems
G2d = ss(-6, 7, 8, 0, 1)
G2d2 = ss(-6, 7, 8, 0, 0.5)
K1d = ss(-1, 1, 1, 0, 1)

# Basic feedback interconnections
@test feedback(K1, ss(1.0)) == ss(-2, 1, 1, 0)
@test feedback(K1, 1.0) == ss(-2, 1, 1, 0)
@test feedback(K1d, ss(1.0, 1)) == ss(-2, 1, 1, 0, 1)
@test_broken feedback(G2d, 1.0) == ss(-2, 1, 1, 0, 1)

# Check that errors for sample-time mismatc are thrown
@test_throws ErrorException feedback(G2, K1d)
@test_throws ErrorException feedback(G2d2, K1d) # Sample time mismatch

# Test general feedback interconnections
@test feedback(G1, K1, U1=[1], Y1=[1], W1=[2], Z1=[2]) == ss([-9 -2; 4 -1], [3; 0], [5 0], 0)
@test feedback(G1, K1, U1=[1], Y1=[1], W1=[1], Z1=[1]) == ss([-9 -2; 4 -1], [2; 0], [4 0], 0)
@test feedback(G1, K1, U1=[2], Y1=[2], W1=[1], Z1=[1]) == ss([-9 -3; 5 -1], [2; 0], [4 0], 0)

@test feedback(K1, G1, W1=Int[], Z1=Int[], U2=[1], Y2=[1], W2=[2], Z2=[2]) == ss([-1 -4; 2 -9], [0; 3], [0 5], 0)

# Feedback with scalar
@test feedback(G2, 1, pos_feedback=false) == ss(-62, 7, 8, 0)
@test feedback(G2, 1, pos_feedback=true) == ss(50, 7, 8, 0)
@test feedback(1, G2, pos_feedback=false) == ss(-62, 7, -8, 1)
@test feedback(1, G2, pos_feedback=true) == ss(50, 7, 8, 1)

# Feedback with scalar (including direct term)
@test feedback(G3, 0.5, pos_feedback=false) ≈ ss(-4/3, 2/3, 2/3, 2/3)
@test feedback(G3, 0.5, pos_feedback=true) == ss(0, 2, 2, 2)
@test feedback(0.5, G3, pos_feedback=false) ≈ ss(-4/3, 1/3, -1/3, 1/3)
@test feedback(0.5, G3, pos_feedback=true) ≈ ss(0, 1, 1, 1)

@test_broken feedback(G3, 1) == ss(-1.5, 0.5, 0.5, 0.5) # Old feedback method
@test feedback(G3, 1, pos_feedback=false) == ss(-1.5, 0.5, 0.5, 0.5)

# Test that errors are thrown for mismatched dimensions
@test_throws ErrorException feedback(G1, K1, U1=1:2, Y2=1)
@test_throws ErrorException feedback(G1, K1, Y2=1:2)


# Tests of linear fractional transformations (LFTs)
@test lft(G1, G2) == ss([-9 24; 35 -6], [2; 0], [4 0], 0)
@test lft(G1, G2, :l) == ss([-9 24; 35 -6], [2; 0], [4 0], 0)
@test lft(G1, G2, :u) == ss([-9 16; 28 -6], [3; 0], [5 0], 0)

@test_throws ErrorException lft(G2, G1)
@test_throws ErrorException lft(G2, G1, :l)
@test_throws ErrorException lft(G2, G1, :u)
@test_throws ErrorException lft(G1, G2, :x) # Invalid type of lft

# Redheffer star product
G4 = ss(-6, [7 8], [11; 12], 0)
@test starprod(G1, G4, 1, 1) == ss([-9 33; 35 -6], [2 0; 0 8], [4 0; 0 12], zeros(2,2))



# Feedback2dof

P0 = tf(1.0, [1, 1, 1])
C = pid(kp=1, ki=1, kd=1)
F = tf(1.0, [1,1])
@test feedback2dof(P0, 0*C, F) == P0*F
@test feedback2dof(P0, C, 0*F) == feedback(P0*C)
@test_nowarn feedback2dof(P0, C, F)


end
