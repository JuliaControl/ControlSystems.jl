using MatrixEquations
include("src/LyapTest.jl")

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5

using Random
Random.seed!(0)

bg = BenchmarkGroup()
for T in [Float64, ComplexF64]
    for k=[5,50,500]#[5, 20, 200]
        n = k
        A = randn(T, n, n)
        B = randn(T, n, n)
        C = randn(T, n, n)
        Csym = Matrix(Hermitian(randn(T, n, n)))

        bg["1. sylvc", 1, T, n] = @benchmarkable LyapTest.sylvc($A, $B, $C)
        bg["1. sylvc", 2, T, n] = @benchmarkable MatrixEquations.sylvc($A, $B, $C)
        bg["1. sylvc", 3, T, n] = @benchmarkable sylvester($A, $B, $C)

        bg["3. sylvd", 1, T, n] = @benchmarkable LyapTest.sylvd($A, $B, $C)
        bg["3. sylvd", 2, T, n] = @benchmarkable MatrixEquations.sylvd($A, $B, $C)

        bg["2. lyapc", 1, T, n] = @benchmarkable LyapTest.lyapc($A, $Csym)
        bg["2. lyapc", 2, T, n] = @benchmarkable MatrixEquations.lyapc($A, $Csym)
        bg["2. lyapc", 3, T, n] = @benchmarkable lyap($A, $Csym)

        bg["4. lyapd", 1, T, n] = @benchmarkable LyapTest.lyapd($A, $Csym)
        bg["4. lyapd", 2, T, n] = @benchmarkable MatrixEquations.lyapd($A, $Csym)
    end
end

results = run(bg)


using PrettyTables


sorted_results = sort(collect(results.data), by=x->(x[1][1], x[1][3]==Float64, x[1][2], x[1][4]))

for k=1:4
    data = reshape(sorted_results[9*(k-1) .+ (1:9)], 3, 3)

    println("$(data[1][1][1]), $(data[1][1][3])")
    table_content = Matrix{Any}(fill(NaN, 3, 4))
    table_content[1:3,1] = [x[1][4] for x in sorted_results[1:3]]
    table_content[1:3,2:4] = [string(BenchmarkTools.prettytime(median(x[2]).time), " (", x[2].allocs, " allocs, ", BenchmarkTools.prettymemory(x[2].memory), ")") for x in data]

    pretty_table(table_content, ["N", "LyapTest", "MatrixEquations", "LinearAlgebra"])
end

for k=1:4
    data = reshape(sorted_results[36 + 6*(k-1) .+ (1:6)], 3, 2)

    #println(data)
    println("$(data[1][1][1]), $(data[1][1][3])")
    table_content = Matrix{Any}(fill(NaN, 3, 3))
    table_content[1:3,1] = [x[1][4] for x in sorted_results[1:3]]
    table_content[1:3,2:3] = [string(BenchmarkTools.prettytime(median(x[2]).time), " (", x[2].allocs, " allocs, ", BenchmarkTools.prettymemory(x[2].memory), ")") for x in data]

    pretty_table(table_content, ["N", "LyapTest", "MatrixEquations"])
end


n = 100

for (numtype, T) in [(:real, Float64), (:complex, ComplexF64)]
    A = randn(T, n, n)
    B = randn(T, n, n)
    C = randn(T, n, n)
    Q = Matrix(Hermitian(randn(T, n, n)))

    println("")
    println("sylvc ($numtype)")
    X = LyapTest.sylvc(A, B, C)
    @printf("LyapTest      err=%.2e\n", norm(A*X + X*B - C))
    X = MatrixEquations.sylvc(A, B, C)
    @printf("MatrixEqs     err=%.2e\n", norm(A*X + X*B - C))
    X = LinearAlgebra.sylvester(A, B, C)
    @printf("LinearAlgebra err=%.2e\n", norm(A*X + X*B + C))

    println("\nlyapc ($numtype)")
    X = LyapTest.lyapc(A, Q)
    @printf("LyapTest      err=%.2e\n", norm(A*X + X*A' + Q))
    X = MatrixEquations.lyapc(A, Q)
    @printf("MatrixEqs     err=%.2e\n", norm(A*X + X*A' + Q))
    X = LinearAlgebra.lyap(A, Q)
    @printf("LinearAlgebra err=%.2e\n", norm(A*X + X*A' + Q))

    println("\nsylvd ($numtype)")
    X = LyapTest.sylvd(A, B, C)
    @printf("LyapTest      err=%.2e\n", norm(A*X*B - X - C))
    X = MatrixEquations.sylvd(A, B, C)
    @printf("MatrixEqs     err=%.2e\n", norm(A*X*B + X - C))

    println("\nlyapd ($numtype)")
    X = LyapTest.lyapd(A, Q)
    @printf("LyapTest      err=%.2e\n", norm(A*X*A' - X + Q))
    X = MatrixEquations.lyapd(A, Q)
    @printf("MatrixEqs     err=%.2e\n", norm(A*X*A' - X + Q))
end


using TimerOutputs
const to = TimerOutput()

for (numtype, T) in [(:real, Float64), (:complex, ComplexF64)]
    n = 100
    A = Matrix(schur(randn(T, n, n)).T') # Lower quasi triangular
    B = schur(randn(T, n, n)).T
    C = randn(T, n, n)
    Csym = Matrix(Hermitian(randn(T, n, n)))

    @timeit to string(numtype) begin
    println("_sylvc_schur (:sylv) ($numtype)")

    X_sylvc = LyapTest._sylvc_schur!(A, B, copy(C), Val(:sylv), Val(numtype))
    @printf("  err=%.2e\n", norm(A*X_sylvc + X_sylvc*B - C))

    X_sylvc = copy(C)
    MatrixEquations.sylvcs!(Matrix(A'), B, X_sylvc)
    @printf("  err=%.2e\n", norm(A'*X_sylvc + X_sylvc*B - C))

    println("_sylvd_schur (:sylv) ($numtype)")

    X_sylvd = LyapTest._sylvd_schur!(A, B, copy(C), Val(:sylv), Val(numtype))
    @printf("  err=%.2e\n", norm(A*X_sylvd*B - X_sylvd - C))

    X_sylvd = copy(C)
    MatrixEquations.sylvds!(Matrix(A'), B, X_sylvd)
    @printf("  err=%.2e\n", norm(A'*X_sylvd*B + X_sylvd - C))


    println("_sylvc_schur (:lyap) ($numtype)")

    X_lyapc = LyapTest._sylvc_schur!(A, Matrix(A'), copy(Csym), Val(:lyap), Val(numtype))
    @printf("  err=%.2e\n", norm(A*X_lyapc + X_lyapc*Matrix(A') - Csym))

    X_lyapc = copy(Csym)
    MatrixEquations.lyapcs!(B, X_lyapc)
    @printf("  err=%.2e\n", norm(B*X_lyapc + X_lyapc*Matrix(B') + Csym))

    println("_sylvd_schur (:lyap) ($numtype)")

    X_lyapd = LyapTest._sylvd_schur!(Matrix(B'), B, copy(Csym), Val(:lyap), Val(numtype))
    @printf("  err=%.2e\n", norm(B'*X_lyapd*B - X_lyapd - Csym))

    X_lyapd = copy(Csym)
    MatrixEquations.lyapds!(B, X_lyapd)
    @printf("  err=%.2e\n", norm(B*X_lyapd*B' - X_lyapd + Csym))
    end
end
