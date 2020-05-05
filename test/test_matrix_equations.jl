# Test _schurstructure
R = diagm(0 => [1, 1, 1, 1, 1, 1, 1, 1], -1 => [1, 0, 0, 0, 1, 0, 0])
d0 = [2, 1, 1, 2, 1, 1]
b0 = [1:2, 3:3, 4:4, 5:6, 7:7, 8:8]
@test (d0, b0, 6) == LyapTest._schurstructure(R, Val(:u))
@test (d0, b0, 6) == LyapTest._schurstructure(R', Val(:l))



# Solve Sylvester equations with A = diagm(a) and B = diagm(b)
sylvc_diag = (a, b, C) -> [C[i,j]/(a[i] + b[j]) for i=1:length(a), j=1:length(b)]
sylvd_diag = (a, b, C) -> [C[i,j]/(a[i]*b[j] - 1) for i=1:length(a), j=1:length(b)]


# Basic tests using 1x1 Matrices

a = [2.0]; A = diagm(a)
ac = [2.0 + im]; Ac = diagm(ac)
b = [3.0]; B = diagm(b)
C = fill(5.0, 1, 1)
Cc = fill(complex(5.0), 1, 1)

@test LyapTest._sylvc_schur!(A, B, copy(C), Val(:sylv), Val(:real)) ≈ sylvc_diag(a, b, C)
@test LyapTest._sylvc_schur!(A, B, copy(C), Val(:sylv), Val(:complex)) ≈ sylvc_diag(a, b, C)
@test LyapTest._sylvc_schur!(Ac, B, copy(Cc), Val(:sylv), Val(:complex)) ≈ sylvc_diag(ac, b, Cc)
@test LyapTest._sylvc_schur!(A, A, copy(C), Val(:lyap), Val(:real)) ≈ sylvc_diag(a, a, C)
@test LyapTest._sylvc_schur!(A, A, copy(C), Val(:lyap), Val(:complex)) ≈ sylvc_diag(a, a, C)
@test LyapTest._sylvc_schur!(Ac', Ac, copy(Cc), Val(:lyap), Val(:complex)) ≈ sylvc_diag(conj(ac), ac, Cc)

@test LyapTest._sylvd_schur!(A, B, copy(C), Val(:sylv), Val(:real)) ≈ sylvd_diag(a, b, C)
@test LyapTest._sylvd_schur!(A, B, copy(C), Val(:sylv), Val(:complex)) ≈ sylvd_diag(a, b, C)
@test LyapTest._sylvd_schur!(Ac, B, copy(Cc), Val(:sylv), Val(:complex)) ≈ sylvd_diag(ac, b, Cc)
@test LyapTest._sylvd_schur!(A, copy(A'), copy(C), Val(:lyap), Val(:real)) ≈ sylvd_diag(a, a, C)
@test LyapTest._sylvd_schur!(A, copy(A'), copy(C), Val(:lyap), Val(:complex)) ≈ sylvd_diag(a, a, C)
@test LyapTest._sylvd_schur!(Ac, copy(Ac'), copy(complex(C)), Val(:lyap), Val(:complex)) ≈ sylvd_diag(ac, conj(ac), C)

@test LyapTest.sylvc(A, B, C) ≈ sylvc_diag(a, b, C)
@test LyapTest.sylvc(Ac, B, Cc) ≈ sylvc_diag(ac, b, Cc)
@test LyapTest.lyapc(A, C) ≈ sylvc_diag(a, a, -C)
@test LyapTest.lyapc(Ac, C) ≈ sylvc_diag(ac, conj(ac), -C)

@test LyapTest.sylvd(A, B, C) ≈ sylvd_diag(a, b, C)
@test LyapTest.sylvd(Ac, B, Cc) ≈ sylvd_diag(ac, b, Cc)
@test LyapTest.lyapd(A, C) ≈ sylvd_diag(a, a, -C)
@test LyapTest.lyapd(Ac, C) ≈ sylvd_diag(conj(ac), ac, -C)


# Basic tests with diagonal 2x2 matrices

a = [2.0, 3]; A = diagm(a)
ac = a .+ im; Ac = diagm(ac)
b = [3.0, 5]; B = diagm(b)
C = [1.0 2; 2 1]
Cc = [1.0 2+im; 1+2im 1]
Cherm = [1.0 2+im; 2-im 1]


@test LyapTest._sylvc_schur!(A, B, copy(C), Val(:sylv), Val(:real)) ≈ sylvc_diag(a, b, C)
@test LyapTest._sylvc_schur!(A, B, copy(C), Val(:sylv), Val(:complex)) ≈ sylvc_diag(a, b, C)
@test LyapTest._sylvc_schur!(Ac, B, copy(Cc), Val(:sylv), Val(:complex)) ≈ sylvc_diag(ac, b, Cc)
@test LyapTest._sylvc_schur!(A, A, copy(C), Val(:lyap), Val(:real)) ≈ sylvc_diag(a, a, C)
@test LyapTest._sylvc_schur!(A, A, copy(C), Val(:lyap), Val(:complex)) ≈ sylvc_diag(a, a, C)
@test LyapTest._sylvc_schur!(copy(Ac'), Ac, copy(Cherm), Val(:lyap), Val(:complex)) ≈ sylvc_diag(ac, conj(ac), Cherm)

@test LyapTest._sylvd_schur!(A, B, copy(C), Val(:sylv), Val(:real)) ≈ sylvd_diag(a, b, C)
@test LyapTest._sylvd_schur!(A, B, copy(C), Val(:sylv), Val(:complex)) ≈ sylvd_diag(a, b, C)
@test LyapTest._sylvd_schur!(Ac, B, copy(Cc), Val(:sylv), Val(:complex)) ≈ sylvd_diag(ac, b, Cc)
@test LyapTest._sylvd_schur!(A, copy(A'), copy(C), Val(:lyap), Val(:real)) ≈ sylvd_diag(a, a, C)
@test LyapTest._sylvd_schur!(A, copy(A'), copy(C), Val(:lyap), Val(:complex)) ≈ sylvd_diag(a, a, C)
@test LyapTest._sylvd_schur!(Ac, copy(Ac'), copy(Cherm), Val(:lyap), Val(:complex)) ≈ sylvd_diag(ac, conj(ac), Cherm)

@test LyapTest.sylvc(A, B, C) ≈ sylvc_diag(a, b, C)
@test LyapTest.sylvc(Ac, B, Cc) ≈ sylvc_diag(ac, b, Cc)
@test LyapTest.lyapc(A, C) ≈ sylvc_diag(a, a, -C)
@test LyapTest.lyapc(Ac, Cherm) ≈ sylvc_diag(conj(ac), ac, -Cherm)

@test LyapTest.sylvd(A, B, C) ≈ sylvd_diag(a, b, C)
@test LyapTest.sylvd(Ac, B, Cc) ≈ sylvd_diag(ac, b, Cc)
@test LyapTest.lyapd(A, C) ≈ sylvd_diag(a, conj(a), -C)
@test LyapTest.lyapd(Ac, Cherm) ≈ sylvd_diag(ac, conj(ac), -Cherm)



#
# Further tests with non-diagonal matrices
#

# Helper function for evaluating the C matrix
sylvc_rhs = (A, B, X) -> (A*X + X*B)
sylvd_rhs = (A, B, X) -> (A*X*B - X)
lyapc_rhs = (A, X) -> -(A*X + X*A')
lyapd_rhs = (A, X) -> -(A*X*A' - X)


for (A, B, X) in [(fill(2.0, 1, 1), fill(3.0, 1, 1), fill(1.0, 1, 1)),
                  ([diagm([2.0,3]), diagm([4.0,5]), [1.0 2; 3 4]]),
                  ([-2 0; 2 3], [3.0 4; 0 5], [1.0 2; 3 4])]

    println(size(A))

    Ac = A .* (1 + im) # Keep the structure of A..
    Xsym = X + X'
    Xherm = (X .+ im) + (X .+ im)'

    @test LyapTest._sylvc_schur!(A, B, sylvc_rhs(A, B, X), Val(:sylv), Val(:real)) ≈ X
    @test LyapTest._sylvc_schur!(A, B, sylvc_rhs(A, B, X), Val(:sylv), Val(:complex)) ≈ X
    @test LyapTest._sylvc_schur!(Ac, B, sylvc_rhs(Ac, B, X), Val(:sylv), Val(:complex)) ≈ X
    @test LyapTest._sylvc_schur!(A, copy(A'), -lyapc_rhs(A, Xsym), Val(:lyap), Val(:real)) ≈ Xsym
    @test LyapTest._sylvc_schur!(A, copy(A'), -lyapc_rhs(A, Xsym), Val(:lyap), Val(:complex)) ≈ Xsym
    @test LyapTest._sylvc_schur!(Ac, copy(Ac'), -lyapc_rhs(Ac, Xherm), Val(:lyap), Val(:complex)) ≈ Xherm

    @test LyapTest._sylvd_schur!(A, B, sylvd_rhs(A, B, X), Val(:sylv), Val(:real)) ≈ X
    @test LyapTest._sylvd_schur!(A, B, sylvd_rhs(A, B, X), Val(:sylv), Val(:complex)) ≈ X
    @test LyapTest._sylvd_schur!(Ac, B, sylvd_rhs(Ac, B, X), Val(:sylv), Val(:complex)) ≈ X
    @test LyapTest._sylvd_schur!(A, copy(A'), -lyapd_rhs(A, Xsym), Val(:lyap), Val(:real)) ≈ Xsym
    @test LyapTest._sylvd_schur!(A, copy(A'), -lyapd_rhs(A, Xsym), Val(:lyap), Val(:complex)) ≈ Xsym
    @test LyapTest._sylvd_schur!(Ac, copy(Ac'), -lyapd_rhs(Ac, Xherm), Val(:lyap), Val(:complex)) ≈ Xherm

end
