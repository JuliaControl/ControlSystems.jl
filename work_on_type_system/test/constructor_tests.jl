# SisoRational constructors
SisoRational{Float64}([1, 0], [1, 1, 1])
SisoRational{Float64}([1.0, 0], [1.0, 1, 1])
SisoRational{Int64}([1, 0], [1, 1, 1])
SisoRational{Complex128}([1 + 2im, 0], [1, 1, 1])

SisoRational([1.0 + 2im, 0], [1, 1, 1])

SisoRational(Poly([0, 1.0 + 2im]), Poly([1, 1, 1]))


# SisoZpk constructors
SisoZpk{Float64,Complex128}(complex([1.0, 0, 1]), complex([1.0, 1, 1]), 1.0)
# SisoZpk{Float64,Complex128}([1.0, 0, 1], [1.0, 1, 1], 1.0) Ska det här funka?

SisoZpk([1.0 + 1im, 1 - 1im], [1im, -1im], 1.0 + 2im)
f_zpk = SisoZpk([1.0 + 1im, 1 - 1im], [1im, -1im], 1.0)

# Construct transfer function directly
TransferFunction{SisoRational{Complex128}}(fill(SisoRational([1.0im+1,0],[1.0im,2]),1,1))
TransferFunction{SisoRational{Float64}}(fill(SisoRational([1.0,0],[1.0,2]),1,1))

TransferFunction{SisoZpk{Float64,Complex128}}(fill(f_zpk,1,1))
TransferFunction{SisoZpk{Float64,Complex128}}(fill(SisoZpk([1.0 + 1im, 1 - 1im], [1im, -1im, 2, 3], 1.0),1,1))

# tf
tf([[1, 0, 2]][:,:], [[2, 0, 3]][:,:], 0.1)
tf([[1.0, 0, 2]][:,:], [[2.0, 0, 3]][:,:], 0.1)

@edit tf(1.0,[1,1])

# zpk
zpk([[1.0im, -1im]][:,:], [[-1.0+1im, -1-1im]][:,:], [0.1][:,:])



zero(SisoRational{Int64})
zero(SisoRational{Complex128})
zero(SisoRational{Float64})

one(SisoRational{Int64})
one(SisoRational{Complex128})
one(SisoRational{Float64})

zero(SisoZpk{Int64})
zero(SisoZpk{Complex128})
zero(SisoZpk{Float64})

one(SisoZpk{Int64, Int64})
one(SisoZpk{Float64, Complex128})
one(SisoZpk{Float64, Float64})

# NOTE: Shouldn't be allowed?
one(SisoZpk{Complex128, Float64})

SisoZpk{Int64}(Int64[], Int64[], 1)

SisoRational(Poly([1, 2]), Poly([1, 2, 1]))



# tf should perhaps also accept polynomials
