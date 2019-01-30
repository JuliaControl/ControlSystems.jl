abstract type LTISystem <: AbstractSystem end
+(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::LTISystem) = -(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::LTISystem) = *(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::LTISystem) = /(promote(sys1, sys2)...)
Base.inv(G::LTISystem) = 1/G

"""`issiso(sys)`

Returns `true` if `sys` is SISO, else returns `false`."""
function issiso(sys::LTISystem)
    return ninputs(sys) == 1 && noutputs(sys) == 1
end


"""`iscontinuous(sys)`

Returns `true` if `sys` is continuous, else returns `false`."""
function iscontinuous(sys::LTISystem)
    return sys.Ts == 0
end


"""`isstable(sys)`

Returns `true` if `sys` is stable, else returns `false`."""
function isstable(sys::LTISystem)
    if iscontinuous(sys)
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





function _check_consistent_sampling_time(systems::AbstractVector{LTISystem})
    if !all(s.Ts == Ts for s in systems)
       error("Sampling time mismatch")
   end
end
function _check_consistent_sampling_time(sys1::LTISystem, sys2::LTISystem)
    if sys1.Ts != sys2.Ts
       error("Sampling time mismatch")
   end
end
