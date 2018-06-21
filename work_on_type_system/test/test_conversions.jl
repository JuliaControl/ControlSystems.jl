f_rational_int = SisoRational([1], [1,2])

convert(SisoRational{Float64}, f_rational_int)

G_tf_int = tf([1], [1,2,3])


tf([1], [1,2,3]) + tf([1.0], [1.0,2,3])

promote_rule(typeof(tf([1], [1,2,3])), typeof(tf([1.0], [1.0,2,3])))



sys1 = ss([1], [2], [3], [0])
sys2 = ss(eye(2), [1 0; 0 2], eye(2), zeros(2,2))

typeof(sys1) <: StateSpace{Int64}
typeof(sys2) <: StateSpace{Float64}


promote_rule(StateSpace{Int64}, StateSpace{Float64})
promote_rule(typeof(sys1), typeof(sys2))



convert(TransferFunction{SisoRational{Float64}}, G_tf_int)
