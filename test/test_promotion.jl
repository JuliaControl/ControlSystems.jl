@testset "test_promotion" begin
P = tf(1)
s = ss(1)
TP = typeof(P)
Ts = typeof(s)
@test promote_type(TP,Ts) == Ts

P = tf(1.)
s = ss(1.)
TP = typeof(P)
Ts = typeof(s)
@test promote_type(Ts,TP) == Ts


P = tf(1.)
s = ss(1)
TP = typeof(P)
Ts = typeof(s)
@test promote_type(Ts,TP) == typeof(ss(1.))


P = zpk(1)
s = ss(1)
TP = typeof(P)
Ts = typeof(s)
@test promote_type(Ts,TP) == typeof(ss(1.))

P = zpk(1)
P2 = zpk(1.)
TP = typeof(P)
TP2 = typeof(P2)
@test promote_type(TP,TP2) == TP2


@test promote_type(typeof(zpk(1)), typeof(zpk(1))) == typeof(zpk(1))
@test promote_type(typeof(tf(1)), typeof(tf(1))) == typeof(tf(1))
@test promote_type(typeof(ss(1)), typeof(ss(1))) == typeof(ss(1))
# @test promote_type(typeof(tfg("1")), typeof(tfg("1"))) == typeof(tfg("1"))

@test promote(ss(1), ss(1.)) == (ss(1.), ss(1.))
@test promote(ss(1), tf(1.)) == (ss(1.), ss(1.))
@test promote(ss(1), zpk(1.)) == (ss(1.), ss(1.))
@test promote(zpk(1), tf(1.)) == (zpk(1.), zpk(1.))
@test promote(tf([1 + im]), tf([1.])) == (tf([1. + im]), tf([1. + 0*im]))
@test promote(zpk(1), tf(1)) == (zpk(1.), zpk(1.))
@test promote(ss(1), ss([1.*im])) == (ss([1+0*im]),ss([1.*im]))
@test promote(ss(1), tf([1.*im])) == (ss([1+0*im]),ss([1.*im]))

end
