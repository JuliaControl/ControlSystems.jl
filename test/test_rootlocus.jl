tf2 = tf([1/5,2],[1,1,1])
rlocusplot(tf2)

# https://github.com/JuliaControl/ControlSystemsBase.jl/issues/740
Nroots = -1.0 .+  [-1.732050807568877im, 1.732050807568877im]
Droots =  [0.0, -4.0, -6.0, -0.7 - 0.7141428428542851im, -0.7 + 0.7141428428542851im]
G = zpk(Nroots, Droots, 1.0)
rlocusplot(G, 200)
