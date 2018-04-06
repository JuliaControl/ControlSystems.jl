abstract type LTISystem end
+(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::LTISystem) = -(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::LTISystem) = *(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::LTISystem) = /(promote(sys1, sys2)...)

@doc """`isstable(sys)`

Returns `true` if `sys` is stable, else returns `false`.""" ->
function isstable(sys::LTISystem)
    if sys.Ts == 0
        if any(real.(pole(sys)).>=0)
            return false
        end
    else
        if any(abs.(pole(sys)).>=1)
            return false
        end
    end
    return true
end

@doc """`iscontinuous(sys)`

Returns `true` if `sys` is continuous, else returns `false`.""" ->
function iscontinuous(sys::LTISystem)
    return sys.Ts == 0
end

@doc """`issiso(sys)`

Returns `true` if `sys` is SISO, else returns `false`.""" ->
function issiso(sys::LTISystem)
    return sys.ny == 1 && sys.nu == 1
end
