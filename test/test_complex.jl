@testset "test_complex_systems" begin

C_1 = zpk([], [-1+im], 1.0+1im)
C_2 = zpk([-1+im], [], 1.0+1im)

# Basic arittmetic
@test C_1 + C_1 == zpk([], [-1+im], 2+2im)
@test (1.0+im)*C_1 == zpk([], [-1+im], 2.0im)
@test C_2 + C_2 == zpk([-1+im], [], 2+2im)

@test C_1 * C_1 == zpk([], [-1+im,-1+im], 2.0im)
@test C_2 * C_2 ≈ zpk([-1+im,-1+im], [], 2.0im)

@test im*ss(1) == ss(im)


@test pole(zpk([], [-1+im,-1+im,0], 2.0im)) == [-1+im,-1+im,0]
@test tzero(zpk([-1+im,-1+im,0], [-2], 2.0im)) == [-1+im,-1+im,0]

@test zpk( tf([1.0, 1+im], [1.0, 2+im]) ) == zpk( [-1-im], [-2-im], 1.0+0im)

@test minreal(zpk([-1+im], [-1+im,-1+im],1+0im)) == zpk([], [-1+im],1+0im)
@test minreal(zpk([-1+im, -1+im], [-1+im],1+1im)) == zpk([-1+im], [], 1+1im)


@test_throws AssertionError zpk([-1+im], [-1+im,-1+im],1) #  Given the type of k this should be a real-coefficient system, but poles and zeros don't come in conjugate pairs

@test zpk([-2+im], [-1+im],1+0im)*zpk([], [-1+im],1+0im) == zpk([-2+im], [-1+im, -1+im], 1+0im)
@test zpk([], [-2], 2) + zpk([], [-1], 1) == zpk([-4/3], [-2,-1], 3)

@test tf(zpk([-2+im], [-1+im],1+0im)) == tf([1, 2-im], [1, 1-im])

end
