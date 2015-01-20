# Create a random, stable StateSpace
function rss(nx::Int, nu::Int=1, ny::Int=1, feedthrough::Bool=true)
    Q = randn(nx, nx)
    A = Q*diagm(-100*abs(randn(nx)))*Q'
    B = randn(nx, nu)
    C = randn(ny, nx)
    if feedthrough
        D = randn(ny, nu)
    else
        D = zeros(ny, nu)
    end
    return ss(A, B, C, D)
end

#########################
## State Space Systems ##
#########################

# Naming convention:
# ------------------
# {type}_{dims}_{stable}[_d]
# type: C: Continuous, D: Discrete
# dims: "nxnuny"
# stable: s: stable, m: marginally stable, u: unstable
# feedthrough: append "_d" if `D` is present

# State matrices
a_1_s = [-5]
a_1_u = [5]
a_2_s = [-5 -3; 2 -9]
a_2_m = [3 -9; 4 -3]
a_2_u = [3 -9; 0 -3]
da_1_s = [-0.5]
da_1_u = [-2]
da_2_s = [0.2 -0.8; -0.8 0.07]
da_2_m = [1 0; 0 0.5]
da_2_u = [-2 1; 0 0.5]

# CONTINUOUS
# nx = 1
C_111_s = ss(a_1_s, [2], [3], [0])
C_111_s_d = ss(a_1_s, [2], [3], [4])
C_111_u = ss(a_1_u, [2], [3], [0])
C_111_u_d = ss(a_1_u, [2], [3], [4])

# nx = 2
C_211_s = ss(a_2_s, [1; 2], [1 0], [0])
C_212_s = ss(a_2_s, [1; 2], eye(2), [0; 0])
C_221_s = ss(a_2_s, [1 0; 0 2], [1 0], [0 0])
C_222_s = ss(a_2_s, [1 0; 0 2], eye(2), zeros(2,2))
C_222_s_d = ss(a_2_s, [1 0; 0 2], eye(2), eye(2))

C_211_m = ss(a_2_m, [1; 2], [1 0], [0])
C_212_m = ss(a_2_m, [1; 2], eye(2), [0; 0])
C_221_m = ss(a_2_m, [1 0; 0 2], [1 0], [0 0])
C_222_m = ss(a_2_m, [1 0; 0 2], eye(2), zeros(2,2))
C_222_m_d = ss(a_2_m, [1 0; 0 2], eye(2), eye(2))

C_222_u = ss(a_2_u, [1 0; 0 2], eye(2), zeros(2,2))

# Pure Gain
C_011 = ss([4])
C_022 = ss(4*eye(2))

# DISCRETE
# nx = 1
D_111_s = ss(da_1_s, [2], [3], [0], 0.005)
D_111_s_d = ss(da_1_s, [2], [3], [4], 0.005)
D_111_u = ss(da_1_u, [2], [3], [0], 0.005)
D_111_u_d = ss(da_1_u, [2], [3], [4], 0.005)

# nx = 2
D_211_s = ss(da_2_s, [1; 2], [1 0], [0], 0.005)
D_212_s = ss(da_2_s, [1; 2], eye(2), [0; 0], 0.005)
D_221_s = ss(da_2_s, [1 0; 0 2], [1 0], [0 0], 0.005)
D_222_s = ss(da_2_s, [1 0; 0 2], eye(2), zeros(2,2), 0.005)
D_222_s_d = ss(da_2_s, [1 0; 0 2], eye(2), eye(2), 0.005)

D_211_m = ss(da_2_m, [1; 2], [1 0], [0], 0.005)
D_212_m = ss(da_2_m, [1; 2], eye(2), [0; 0], 0.005)
D_221_m = ss(da_2_m, [1 0; 0 2], [1 0], [0 0], 0.005)
D_222_m = ss(da_2_m, [1 0; 0 2], eye(2), zeros(2,2), 0.005)
D_222_m_d = ss(da_2_m, [1 0; 0 2], eye(2), eye(2), 0.005)

D_222_u = ss(da_2_u, [1 0; 0 2], eye(2), zeros(2,2), 0.005)

# Pure Gain
D_011 = ss([4], 0.005)
D_022 = ss(4*eye(2), 0.005)
