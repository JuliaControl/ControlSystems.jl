@testset "test_promotion" begin
G = tf(1)
sys = ss(1)
TG = typeof(G)
Tsys = typeof(sys)
@test promote_type(Tsys,TG) == Tsys

G = tf(1.)
sys = ss(1.)
TG = typeof(G)
Tsys = typeof(sys)
@test promote_type(Tsys,TG) == Tsys


G = tf(1.)
sys = ss(1)
TG = typeof(G)
Tsys = typeof(sys)
@test promote_type(Tsys,TG) == typeof(ss(1.))


G = zpk(1)
sys = ss(1)
TG = typeof(G)
Tsys = typeof(sys)
@test promote_type(Tsys,TG) == typeof(ss(1))

G1 = zpk(1)
G2 = zpk(1.)
TG1 = typeof(G1)
TG2 = typeof(G2)
@test promote_type(TG1,TG2) == TG2

@test promote_type(typeof(zpk(1)), typeof(zpk(1))) == typeof(zpk(1))
@test promote_type(typeof(tf(1)), typeof(tf(1))) == typeof(tf(1))
@test promote_type(typeof(ss(1)), typeof(ss(1))) == typeof(ss(1))


@test promote(ss(1), 1.0) == (ss(1.), ss(1.))
@test promote(ss(1), 1) == (ss(1.), ss(1.))
@test promote(ss(1), ss(1.)) == (ss(1.), ss(1.))
@test promote(ss(1), tf(1.)) == (ss(1.), ss(1.))
@test promote(ss(1), zpk(1.)) == (ss(1.), ss(1.))
@test promote(zpk(1), tf(1.)) == (zpk(1.), zpk(1.))
@test promote(tf([1 + im]), tf([1.])) == (tf([1. + im]), tf([1. + 0*im]))
@test promote(zpk(1), tf(1)) == (zpk(1.), zpk(1.))
@test promote(ss(1), ss([1 .*im])) == (ss([1+0*im]),ss([1 .*im]))
@test promote(ss(1), tf([1 .*im])) == (ss([1+0*im]),ss([1 .*im]))
@test promote(tf(1), 1) == (tf(1.), tf(1.))


end
