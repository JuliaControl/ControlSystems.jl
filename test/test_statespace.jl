
# Naming convention:
# ------------------
# {type}[S]_{dims}[_d][_n]
# type: C: Continuous, D: Discrete
# scalar: add "S" if any input is scalar
# dims: "nxnuny"
# feedthrough: append "_d" if `D` is present
# names: append "_n" if some inputs/outputs/states are named
@testset "test_statespace" begin
    # SCALARS
    for SS in (StateSpace, HeteroStateSpace), SS2 in (StateSpace, HeteroStateSpace)
        a_2 = [-5 -3; 2 -9]
        CS_111 = @inferred SS(-5, 2, 3, [0])
        CS_111_d = @inferred SS([3], 2, 1, 1)
        CS_211 = @inferred SS(a_2, [1; 2], [1 0], 0)
        CS_221 = @inferred SS(a_2, [1 0; 0 2], [1 0], 0)
        CS_222 = @inferred SS(a_2, [1 0; 0 2], eye_(2), 0)

        # CONTINUOUS
        a_1 = [-5]
        C_111 = @inferred SS(a_1, [2], [3], [0])
        C_211 = @inferred SS2(a_2, [1; 2], [1 0], [0])
        C_212 = @inferred SS2(a_2, [1; 2], eye_(2), [0; 0])
        C_221 = @inferred SS2(a_2, [1 0; 0 2], [1 0], [0 0])
        C_222 = @inferred SS(a_2, [1 0; 0 2], eye_(2), zeros(Int,2,2))
        C_222_d = @inferred SS(a_2, [1 0; 0 2], eye_(2), eye_(2))
        C_022 = @inferred SS(4.0*eye_(2))

        # DISCRETE
        da_1 = [-0.5]
        da_2 = [0.2 -0.8; -0.8 0.07]
        D_111 = @inferred SS(da_1, [2], [3], [0], 0.005)
        D_211 = @inferred SS2(da_2, [1; 2], [1 0], [0], 0.005)
        D_221 = @inferred SS2(da_2, [1 0; 0 2], [1 0], [0 0], 0.005)
        D_222 = @inferred SS(da_2, [1 0; 0 2], eye_(2), zeros(2,2), 0.005)
        D_222_d = @inferred SS(da_2, [1 0; 0 2], eye_(2), eye_(2), 0.005)
        D_022 = @inferred SS(4.0*eye_(2), 0.005)

        # Definition of input, output and state names
        C_222_d_n = SS(a_2, [1 0; 0 2], eye_(2), eye_(2))

        # TESTS
        # Contstuct with scalars
        @test CS_111 == C_111
        @test CS_111_d == SS([3],[2],[1],[1])
        @test CS_211 == C_211
        @test CS_221 == C_221
        @test CS_222 == C_222

        # Addition
        @test C_111 + C_111 == SS([-5 0; 0 -5],[2; 2],[3 3],[0])
        @test C_222 + C_222 == SS([-5 -3 0 0; 2 -9 0 0; 0 0 -5 -3;
        0 0 2 -9],[1 0; 0 2; 1 0; 0 2], [1 0 1 0; 0 1 0 1],[0 0; 0 0])
        @test C_222 + 1 == SS([-5 -3; 2 -9],[1 0; 0 2],[1 0; 0 1],[1 1; 1 1])
        @test D_111 + D_111 == SS([-0.5 0; 0 -0.5],[2; 2],[3 3],[0], 0.005)

        @inferred C_111 + C_111

        @test C_111 + false == C_111
        @test false + C_111 == C_111
        @test 1.0*C_111 + false == C_111

        @inferred C_111 + false

        @test C_222 + 1.5 == 1.0C_222 + 1.5 # C_222 has eltype Int
        @test 1.5 + C_222 == 1.0C_222 + 1.5

        # Subtraction
        @test C_111 - C_211 == SS([-5 0 0; 0 -5 -3; 0 2 -9],[2; 1; 2],[3 -1 -0],[0])
        @test 1 - C_222 == SS([-5 -3; 2 -9],[1 0; 0 2],[-1 -0; -0 -1],[1 1; 1 1])
        @test D_111 - D_211 == SS([-0.5 0 0; 0 0.2 -0.8; 0 -0.8 0.07],[2; 1; 2],
        [3 -1 -0],[0], 0.005)

        # Multiplication
        @test C_111 * C_221 == SS([-5 2 0; 0 -5 -3; 0 2 -9],
        [0 0; 1 0; 0 2],[3 0 0],[0 0])
        @test C_212 * C_111 == SS([-5 -3 3; 2 -9 6; 0 0 -5],
        [0; 0; 2],[1 0 0; 0 1 0],[0; 0])
        @test 4*C_222 == SS([-5 -3; 2 -9],[1 0; 0 2],[4 0; 0 4],[0 0; 0 0])
        @test D_111 * D_221 == SS([-0.5 2 0; 0 0.2 -0.8; 0 -0.8 0.07],
        [0 0; 1 0; 0 2],[3 0 0],[0 0],0.005)
        @test C_111 * I(2) == I(2) * C_111 == SS(diagm([a_1; a_1]), 2*I(2), 3*I(2), 0*I(2))
        @test minreal(C_111*C_222_d - C_222_d*C_111, atol=1e-3) == ss(0*I(2)) # scalar times MIMO
        @test C_111*C_222 == ss([-5 0 2 0; 0 -5 0 2; 0 0 -5 -3; 0 0 2 -9], [0 0; 0 0; 1 0; 0 2], [3 0 0 0; 0 3 0 0], 0)

        @inferred C_111 * C_221
        @inferred C_111 * I(2)

        # Test that multiplication/division is applied at correct input/output location
        @test (10*C_111).C == 10*C_111.C
        @test (10*C_111).B == C_111.B

        @test (C_111*10).C == C_111.C
        @test (C_111*10).B == 10*C_111.B

        @test (10\C_111).C == 10\C_111.C
        @test (10\C_111).B == C_111.B

        @test (C_111/10).C == C_111.C
        @test (C_111/10).B == C_111.B/10

        # Division
        @test 1/C_222_d == SS([-6 -3; 2 -11],[1 0; 0 2],[-1 0; -0 -1],[1 -0; 0 1])
        @test C_221/C_222_d == SS([-5 -3 -1 0; 2 -9 -0 -2; 0 0 -6 -3;
        0 0 2 -11],[1 0; 0 2; 1 0; 0 2],[1 0 0 0],[0 0])
        @test 1/D_222_d == SS([-0.8 -0.8; -0.8 -1.93],[1 0; 0 2],[-1 0; -0 -1],
        [1 -0; 0 1],0.005)

        fsys = ss(1,1,1,0)/3 # Int becomes FLoat after division
        @test fsys.B[]*fsys.C[] == 1/3

        @inferred 1/C_222_d
        @inferred C_221/C_222_d

        # Indexing
        @test size(C_222) == (2, 2)
        @test size(C_212) == (2, 1)
        @test C_222[1,1] == SS([-5 -3; 2 -9],[1; 0],[1 0],[0])
        @test C_222[1:1,1] == SS([-5 -3; 2 -9],[1; 0],[1 0],[0])
        @test C_222[1,1:2] == C_221
        @test size(C_222[1,[]]) == (1,0)
        @test C_222[end, end] == C_222[2,2]


        A = [-1.0 -2.0; 0.0 -1.0]
        B = [0.0; -2.0]
        C = [1.0 1.0]
        D = 1.0
        sys = SS(A, B, C, D)

        @test sys + 1.0 == SS(A, B, C, D + 1.0)
        @test 2.0 + sys == SS(A, B, C, D + 2.0)

        @test -sys == SS(A, B, -C, -D)

        ## Test specifying time using TimeEvolution struct
        @test SS(A, B, C, D) == SS(A, B, C, D, Continuous())
        @test SS(A, B, C, D, 1.5) == SS(A, B, C, D, Discrete(1.5))

        # Accessing Ts through .Ts
        @test D_111.Ts == 0.005

        # propertynames
        @test propertynames(C_111) == (:A, :B, :C, :D, :timeevol, :nu, :ny, :nx)
        @test propertynames(D_111) == (:A, :B, :C, :D, :timeevol, :nu, :ny, :nx, :Ts)

        # Printing
        if SS <: StateSpace
            @test sprint(show, C_222) == "StateSpace{Continuous, Int64}\nA = \n -5  -3\n  2  -9\nB = \n 1  0\n 0  2\nC = \n 1  0\n 0  1\nD = \n 0  0\n 0  0\n\nContinuous-time state-space model"
            @test sprint(show, C_022) == "StateSpace{Continuous, Float64}\nD = \n 4.0  0.0\n 0.0  4.0\n\nContinuous-time state-space model"
            @test sprint(show, D_022) == "StateSpace{Discrete{Float64}, Float64}\nD = \n 4.0  0.0\n 0.0  4.0\n\nSample Time: 0.005 (seconds)\nDiscrete-time state-space model"
            @test sprint(show, D_222) == "StateSpace{Discrete{Float64}, Float64}\nA = \n  0.2  -0.8\n -0.8   0.07\nB = \n 1.0  0.0\n 0.0  2.0\nC = \n 1.0  0.0\n 0.0  1.0\nD = \n 0.0  0.0\n 0.0  0.0\n\nSample Time: 0.005 (seconds)\nDiscrete-time state-space model"
        end

        # Different types
        if VERSION >= v"1.8.0-rc1"
            K1 = SS(I(2)) # Bool
            K2 = SS(1.0I(2)) # Float64
            P = ssrand(3,3,2)
            @test lft(P, -K1) == lft(P, -K2)
        else # With HeteroStateSpace, we get an error on julia < 1.8 due to inverse of sparse I
            K1 = ss(I(2)) # Bool
            K2 = ss(1.0I(2)) # Float64
            P = ssrand(3,3,2)
            @test lft(P, -K1) == lft(P, -K2)
        end

        # Errors
        @test_throws ErrorException C_111 + C_222             # Dimension mismatch
        @test_throws ErrorException C_111 - C_222             # Dimension mismatch
        @test_throws ErrorException D_111 + C_111             # Sampling time mismatch
        @test_throws ErrorException D_111 - C_111             # Sampling time mismatch
        @test_throws ErrorException D_111 * C_111             # Sampling time mismatch
        D_diffTs = SS([1], [2], [3], [4], 0.1)
        @test_throws ErrorException D_111 + D_diffTs            # Sampling time mismatch
        @test_throws ErrorException D_111 - D_diffTs            # Sampling time mismatch
        @test_throws ErrorException D_111 * D_diffTs            # Sampling time mismatch
        @test_throws ErrorException 1/C_222                     # Not invertible
        @test_throws ErrorException 1/C_212                     # Not invertible
        @test_throws ErrorException SS([1 2], [1], [2], [3])      # Not square A
        @test_throws ErrorException SS([1], [2 0], [1], [2])      # I/0 dim mismatch
        @test_throws ErrorException SS([1], [2], [3 4], [1])      # I/0 dim mismatch
        @test_throws ErrorException SS([1], [2], [3], [4], -0.1)  # Negative samping time
        @test_throws ErrorException SS(eye_(2), eye_(2), eye_(2), [0]) # Dimension mismatch

        # Misc tests
        @test SS(1) == SS(zeros(Int,0,0),zeros(Int,0,1),zeros(Int, 1,0),1)
        @test SS(1, 0.1) == SS(zeros(Int,0,0),zeros(Int,0,1),zeros(Int, 1,0),1, 0.1)
        @test zero(StateSpace{Continuous, Float64}) == SS(0)
        @test ss(1) - 1 == (1 - ss(1)) == zero(StateSpace{Continuous, Int})
        @test ndims(C_111) == 2
        @test_nowarn print(C_111)
        @test isproper(C_222)
        @test !isproper(C_222_d)

        @test SS(1, 0.1).Ts == 0.1
        if VERSION >= v"1.8.0-rc1"
            @test @test_logs (:warn, r"deprecated") SS(1).Ts == 0
        end

    end

end
