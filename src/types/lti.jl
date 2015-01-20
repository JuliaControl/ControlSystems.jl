abstract LTISystem

+(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
-(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
*(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
/(sys1::LTISystem, sys2::LTISystem) = +(promote(sys1, sys2)...)
