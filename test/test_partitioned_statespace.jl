using ControlSystems
using LinearAlgebra

##
# system of the form
# 1 | 2 3
# ———————
# 4 | 5 6
# 7 | 8 9
sys1 = ControlSystems.PartionedStateSpace(ss(1.0, [2.0 3.0], [4.0; 7.0], [5.0 6.0; 8.0 9.0]), 1, 1)

matrix(x::Number) = fill(x, 1, 1)

sys1.A == matrix(1.0)
sys1.B1 == matrix(2.0)
sys1.B2 == matrix(3.0)
sys1.C1 == matrix(4.0)
sys1.C2 == matrix(7.0)
sys1.D11 == matrix(5.0)
sys1.D12 == matrix(6.0)
sys1.D21 == matrix(8.0)
sys1.D22 == matrix(9.0)

##
sys2 = ControlSystems.PartionedStateSpace(ss(fill(1.0, 2, 2), fill(2.0, 2, 5), fill(3.0, 7, 2), fill(4.0, 7, 5)), 2, 3)

@test sys2.A == fill(1.0, 2, 2)
@test sys2.B1 == fill(2.0, 2, 2)
@test sys2.B2 == fill(2.0, 2, 3)
@test sys2.C1 == fill(3.0, 3, 2)
@test sys2.C2 == fill(3.0, 4, 2)
@test sys2.D11 == fill(4.0, 3, 2)
@test sys2.D12 == fill(4.0, 3, 3)
@test sys2.D21 == fill(4.0, 4, 2)
@test sys2.D22 == fill(4.0, 4, 3)



@test sys1 == deepcopy(sys1)
@test sys1 != deepcopy(sys2)

# TODO: Add some tests for interconnections, implicitly tested through delay system implementations though
@test (sys1 + sys1).P[1, 1] == (sys1.P[1,1] + sys1.P[1,1])

@test propertynames(sys1) == (:P, :nu1, :ny1)
