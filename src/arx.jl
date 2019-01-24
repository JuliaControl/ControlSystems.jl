"""
    getARXregressor(y::AbstractVector,u::AbstractVecOrMat, na, nb)
Returns a shortened output signal `y` and a regressor matrix `A` such that the least-squares ARX model estimate of order `na,nb` is `y\\A`
Return a regressor matrix used to fit an ARX model on, e.g., the form
`A(z)y = B(z)f(u)`
with output `y` and input `u` where the order of autoregression is `na` and
the order of input moving average is `nb`
# Example
Here we test the model with the Function `f(u) = √(|u|)`
```julia
A     = [1,2*0.7*1,1] # A(z) coeffs
B     = [10,5] # B(z) coeffs
u     = randn(100) # Simulate 100 time steps with Gaussian input
y     = filt(B,A,u)
yr,A  = getARXregressor(y,u,3,2) # We assume that we know the system order 3,2
x     = A\\yr # Estimate model polynomials
plot([yr A*x], lab=["Signal" "Prediction"])
```
For nonlinear ARX-models, see [BasisFunctionExpansions.jl](https://github.com/baggepinnen/BasisFunctionExpansions.jl/)
"""
function getARXregressor(y::AbstractVector,u::AbstractVecOrMat, na, nb)
    length(nb) == size(u,2) || throw(ArgumentError("Length of nb must equal number of input signals"))
    m    = max(na+1,maximum(nb)) # Start of yr
    n    = length(y) - m + 1 # Final length of yr
    A    = toeplitz(y[m:m+n-1],y[m:-1:1])
    y    = A[:,1] # extract yr
    A    = A[:,2:end]
    for i = 1:length(nb)
        s = m - nb[i]
        offs = s - nb[i]+1
        A = [A toeplitz(u[s:s+n-1,i],u[s:-1:offs,i])]
    end
    return y,A
end

"""
    find_na(y::AbstractVector,n::Int)
Plots the RMSE and AIC For model orders up to `n`. Useful for model selection
"""
function find_na(y::AbstractVector,n::Int)
    error = zeros(n,2)
    for i = 1:n
        w,e = ar(y,i)
        error[i,1] = rms(e)
        error[i,2] = aic(e,i)
        print(i,", ")
    end
    println("Done")
    scatter(error, show=true)
end

rms(x::AbstractVector) = sqrt(mean(abs2,x))
sse(x::AbstractVector) = x⋅x

rms(x::AbstractMatrix) = sqrt.(mean(abs2.(x),dims=2))[:]
sse(x::AbstractMatrix) = sum(abs2,x,dims=2)[:]
modelfit(y,yh) = 100 * (1 .-rms(y.-yh)./rms(y.-mean(y)))
aic(x::AbstractVector,d) = log(sse(x)) .+ 2d/size(x,2)
const nrmse = modelfit

"""
    Gtf, Σ = arx(h,y, u, na, nb; λ = 0)

Fit a transfer Function to data using an ARX model and equation error minimization.
`nb` and `na` are the length of the numerator and denominator polynomials. `h` is the sample time of the data. `λ > 0` can be provided for L₂ regularization.
The number of free parameters is `na-1+nb`
`Σ` is the covariance matrix of the parameter estimate. See `bodeconfidence` for visualiztion of uncertainty.

Supports MISO estimation by supplying a matrix `u` where times is first dim, with nb = [nb₁, nb₂...]
"""
function arx(h,y::AbstractVector, u::AbstractVecOrMat, na, nb; λ = 0)
    all(nb .< na) || throw(DomainError(nb,"nb must be smaller than na"))
    na >= 1 || throw(ArgumentError("na must be positive"))
    na -= 1
    y_train, A = getARXregressor(y,u, na, nb)

    if λ == 0
        w = A\y_train
    else
        w = (A'A + λ*I)\A'y_train
    end
    a,b = params2poly(w,na,nb)
    model = tf(b,a,h)
    Σ = parameter_covariance(y_train, A, w, λ)
    return model, Σ
end

# Helper constructor to make a MISO system after MISO arx estimation
function tf(b::AbstractVector{<:AbstractVector{<:Number}}, a::AbstractVector{<:Number}, h)
    tfs = map(b) do b
        tf(b,a,h)
    end
    hcat(tfs...)
end

"""
    a,b = params2poly(params,na,nb)
Used to get numerator and denominator polynomials after arx fitting
"""
function params2poly(w,na,nb)
    a = [1; -w[1:na]]
    w = w[na+1:end]
    b = map(nb) do nb
        b = w[1:nb]
        w = w[nb+1:end]
        b
    end
    a,b
end

"""
    Σ = parameter_covariance(y_train, A, w, λ=0)
"""
function parameter_covariance(y_train, A, w, λ=0)
    σ² = var(y_train .- A*w)
    iATA = if λ == 0
        inv(A'A)
    else
        ATA = A'A
        ATAλ = factorize(ATA + λ*I)
        ATAλ\ATA/ATAλ
    end
    iATA = (iATA+iATA')/2
    Σ = σ²*iATA + sqrt(eps())*Matrix(LinearAlgebra.I,size(iATA))
end

"""
    bodeconfidence(arxtf::TransferFunction, Σ::Matrix; ω = logspace(0,3,200))
Plot a bode diagram of a transfer function estimated with [`arx`](@ref) with confidence bounds on magnitude and phase.
"""
bodeconfidence

@userplot BodeConfidence

@recipe function BodeConfidence(p::BodeConfidence; ω = exp10.(LinRange(-2,3,200)))
    arxtfm = p.args[1]
    Σ      = p.args[2]
    L      = cholesky(Hermitian(Σ)).L
    am, bm = -denpoly(arxtfm)[1].a[2:end], arxtfm.matrix[1].num.a
    wm     = [am; bm]
    na,nb  = length(am), length(bm)
    mc     = 100
    res = map(1:mc) do _
        w             = L*randn(size(L,1)) .+ wm
        a,b           = params2poly(w,na,nb)
        arxtf         = tf(b,a,arxtfm.Ts)
        mag, phase, _ = bode(arxtf, ω)
        mag[:], phase[:]
    end
    magmc      = reduce(hcat, getindex.(res,1))
    phasemc    = reduce(hcat, getindex.(res,2))
    mag        = mean(magmc,dims=2)[:]
    phase      = mean(phasemc,dims=2)[:]
    # mag,phase,_ = bode(arxtfm, ω) .|> x->x[:]
    uppermag   = getpercentile(magmc,0.95)[:]
    lowermag   = getpercentile(magmc,0.05)[:]
    upperphase = getpercentile(phasemc,0.95)[:]
    lowerphase = getpercentile(phasemc,0.05)[:]
    layout := (2,1)

    @series begin
        subplot := 1
        title --> "ARX estimate"
        ylabel --> "Magnitude"
        fillrange := (lowermag, uppermag)
        yscale --> :log10
        xscale --> :log10
        alpha --> 0.3
        ω, mag
    end
    @series begin
        subplot := 2
        fillrange := (lowerphase, upperphase)
        ylabel --> "Phase [deg]"
        xlabel --> "Frequency [rad/s]"
        xscale --> :log10
        alpha --> 0.3
        ω, phase
    end

end

"""
    getpercentile(x,p)

calculates the `p`th percentile along dim 2
"""
function getpercentile(mag,p)
    uppermag = mapslices(mag, dims=2) do magω
        sort(magω)[round(Int,length(magω)*p)]
    end
end
