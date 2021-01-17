@testset "test_autovec" begin


# Test all different autovecd methods
sys = tf([1], [1, 2])
w = exp10.(LinRange(-3,3,10))

# Check output of bode and make sure dimensions are correct
mag, phase, w = bode(sys)
@test size(mag) == size(phase) == (size(w,1), 1, 1)

# Check output of bodev and make sure dimensions are correct
magv, phasev, w = bodev(sys)
@test size(magv) == size(phasev) == size(w)
@test vec(mag) == magv && vec(phase) == phasev

mag, phase, _ = bodev(sys, w)
@test size(mag) == size(phase) == size(w)

# Check output of nyquist and make sure dimensions are correct
real, imag, w = nyquist(sys)
@test size(real) == size(imag) == (size(w,1), 1, 1)

# Check output of nyquistv and make sure dimensions are correct
realv, imagv, w = nyquistv(sys)
@test size(realv) == size(imagv) == size(w)
@test vec(real) == realv && vec(imag) == imagv

real, imag, _ = nyquistv(sys, w)
@test size(real) == size(imag) == size(w)

# Check output of sigma and make sure dimensions are correct
sv, w = sigma(sys)
@test size(sv) == (size(w,1), 1)

# Check output of sigmav and make sure dimensions are correct
svv, w = sigmav(sys)
@test size(svv) == size(w)
@test vec(sv) == svv 

sv, _ = sigmav(sys, w)
@test size(sv) == size(w)

# Check output of freqresp and make sure dimensions are correct
fs = freqresp(sys, w)
@test size(fs) == (size(w,1), 1, 1)

# Check output of freqrespv and make sure dimensions are correct
fsv = freqrespv(sys, w)
@test size(fsv) == size(w)
@test vec(fs) == fsv 






# Test that we can define varous kinds of methods with autovec

# Test arguments for single output
ControlSystems.@autovec () function t0(b::Int, c::Float64=1.0)
    return [b c]
end

@test @isdefined t0
@test @isdefined t0v
@test t0(5) == [5.0 1.0]
@test t0v(5) == [5.0, 1.0]

# Test arguments
ControlSystems.@autovec (2,) function t0(a, b::Int, c::Float64=1.0)
    return a, [b c]
end

@test @isdefined t0
@test @isdefined t0v
@test t0(4, 5) == (4,[5.0 1.0])
@test t0v(4, 5) == (4,[5.0, 1.0])

@test t0v(4, 5, 6.0) == (4,[5.0, 6.0])
@test t0v(4, 5, 6.0) == (4,[5.0, 6.0])

# Test keyword arguments
ControlSystems.@autovec (2,) function t1(; kw=1)
    return kw, [kw kw]
end

@test @isdefined t1
@test @isdefined t1v
@test t1() == (1,[1 1])
@test t1v() == (1,[1, 1])

@test t1(kw=2) == (2,[2 2])
@test t1v(kw=2) == (2,[2, 2])

# Test type constraints on keyword arguments
ControlSystems.@autovec (2,) function t2(; kw::T=1) where {T <: Integer}
    return kw, [kw kw]
end

@test @isdefined t2
@test @isdefined t2v
@test t2() == (1,[1 1])
@test t2v() == (1,[1, 1])

@test t2(kw=2) == (2,[2 2])
@test t2v(kw=2) == (2,[2, 2])

@test_throws MethodError t2(kw=2.) # This method should not be defined

# Test arguments simple method definition
ControlSystems.@autovec (2,) t3(a, b::Int, c::Float64=1.0) = a, [b c]

@test @isdefined t3
@test @isdefined t3v
@test t3(4, 5) == (4,[5.0 1.0])
@test t3v(4, 5) == (4,[5.0, 1.0])

@test t3v(4, 5, 6.0) == (4,[5.0, 6.0])
@test t3v(4, 5, 6.0) == (4,[5.0, 6.0])

# Test keyword arguments simple method definition
ControlSystems.@autovec (2,) t4(; kw=1) = kw, [kw kw]

@test @isdefined t4
@test @isdefined t4v
@test t4() == (1,[1 1])
@test t4v() == (1,[1, 1])

@test t4(kw=2) == (2,[2 2])
@test t4v(kw=2) == (2,[2, 2])

# Test type constraints on keyword arguments simple method definition
ControlSystems.@autovec (2,) t5(; kw::T=1) where {T <: Integer} = kw, [kw kw]

@test @isdefined t5
@test @isdefined t5v
@test t5() == (1,[1 1])
@test t5v() == (1,[1, 1])

@test t5(kw=2) == (2,[2 2])
@test t5v(kw=2) == (2,[2, 2])

@test_throws MethodError t2(kw=2.) # This method should not be defined

# Check MIMO
mimo_sys = ss([-1 1; 0 -3], [0 1; 1 2], [0 1], [0 0])
@test_throws ArgumentError bodev(mimo_sys) # This method should throw error



end
