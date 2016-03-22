using ControlSystems

# function time2freqdomain(ω, lowerbound, upperbound, timevec=linspace(0,10,100), w = [0.1:0.25:4.5; 5:1:10], zeta = [0:0.025:1; 1.2:0.2:2])
function time2freqdomain(ω, lowerbound, upperbound, timevec=linspace(0,10,100), w = logspace(-3,3,30), zeta = logspace(-2,1,30))
  fig1 = Plots.plot()
  # fig2 = Plots.plot()
  ampmax = Float64[]
  ampmin = Float64[]
  for w = w
    for zeta = zeta
      g = tf(w^2,[1,2w*zeta, w^2])
      y,t = step(g, timevec)[1:2]
      if maximum(y-upperbound) <=0 && minimum(y-lowerbound) >= 0
        Plots.plot!(fig1,t,y, c=:blue)
        Plots.gui()
        amp = abs(squeeze(freqresp(g,2*pi*ω)[1],tuple(1,2)))
        ampmax = maximum([ampmax;amp])
        ampmin = minimum([ampmin;amp])
        # Plots.plot!(fig2, ω,abs(amp), c=:blue, xscale=:log10,  yscale=:log10)
      end
    end
  end
  Plots.plot!(fig1,timevec,[lowerbound upperbound], c=:black)
end

timevec = linspace(0,10,1000)
lowerbound = [linspace(-0.1,0.7,300); 0.7*ones(700)]
upperbound = [1.5*ones(200); 1.3*ones(800)]
ω = logspace(-3,3,40)

time2freqdomain(ω, lowerbound, upperbound, timevec)



function calc_sensitivity_bounds()
  fasv = -270:2:0
  for fas = fasv
    for ampdb = -40:0.5:10
      c = 10^((ampdb-ampnom)/20)*exp(im*(fas-fasnom)*pi/180)
      l = template*c
      s = 1./(1+l)
      if maximum(abs(s)) > spec
        bound=[bound ampdb]
        fasbound = [fasbound fas]
        break
      end
    end
  end
  for fas = fasv[end:-1:1]
    for ampdb = 10:-0.5:-20
      c = 10^((ampdb-ampnom)/20)*exp(im*(fas-fasnom)*pi/180)
      l = template*c
      s = 1./(1+l)
      if maximum(abs(s)) > spec
        bound=[bound ampdb]
        fasbound = [fasbound fas]
        break
      end
    end
  end
  bound = [bound bound[1]]
  fasbound = [fasbound fasbound[1]]
  plot(fasbound,bound,c=:red,linewidth=2)
end
