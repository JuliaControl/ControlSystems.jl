v = roots2real_poly_factors([-1.0+1im, -1-1im, 4, 2])


f = SisoZpk{Float64,Complex128}([-1.0 + im,-2,-3], [-5.0 + im, -5.0-im], 1.0)



G = tf(reshape(1:9, 3, 3)[:,:])



tf([1, 2], [1, 2, 1]) + .5 == tf([1, 2.5], [1, 2, 1])



tf([1, 2], [1, 2, 1]) + .5
