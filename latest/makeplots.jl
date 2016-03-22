pwd()
mkdir("docs/build/plots")

Plots.pyplot()

## PID design functions
P = tf(1,[1,1])^4
gangoffourplot(P,tf(1))
Plots.savefig("docs/build/plots/pidgofplot.svg")

## Advanced pole-zero placement
#ζ = 0.2
#ω = 1

#B = [1]
#A   = [1, 2ζ*ω, ω^2]
#P  = tf(B,A)#ζ0 = 0.7
#ω0 = 2
#Am = [1, 2ζ0*ω0, ω0^2]
#Ao = conv(2Am, [1/2, 1])
#AR = [1,0]
#B⁺  = [1]
#B⁻  = [1]
#Bm  = conv(B⁺, B⁻)
#R,S,T = rstc(B⁺,B⁻,A,bm,Am,Ao,AR)
#Gcl = tf(conv(B,T),zpconv(A,R,B,S))
#f1 = stepplot([P,Gcl])
#f2 = gangoffourplot(P, tf(-S,R))
#Plots.savefig(f1, "rstcstepplot.svg")
#Plots.savefig(f2, "rstcgofplot.svg")
