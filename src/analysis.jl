@doc """`pole(sys)`

Compute the poles of system `sys`.""" ->
pole(sys::StateSpace) = eig(sys.A)[1]
pole(sys::TransferFunction) = [map(pole, sys.matrix)...]
