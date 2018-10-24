using Plots
"""
This is just a hack to visualize the plots, should be made into a recepie
eventually and included as a regular plot option.
"""
function visualize_synthesis(Pcl, S, KS, T, gamma; filename=[], tmax=10, fmin=-6, fmax=6)

  # Plot the singular values  of the system
  f = [10^i for i in range(fmin, stop=fmax, length=1001)]

  valPcl  = sigma(Pcl, f)[1];
  valS    = sigma(S, f)[1];
  valKS   = sigma(KS, f)[1];
  valT    = sigma(T, f)[1];

  # Visualize the close loop gain
  pGain = plot(f, valPcl[:,1], xscale = :log10, yscale = :log10, color = :black, w=1, label="\$\\sigma(P_{w\\rightarrow z}(j\\omega))\$")
  plot!(f, gamma*ones(size(f)), xscale = :log10, yscale = :log10, color = :black, w=2, style=:dot, label="\$\\gamma\$")
  if size(valPcl,2) > 1
    for ii = 2:size(valPcl,2)
      plot!(f, valPcl[:,ii], xscale = :log10, yscale = :log10, color = :black, w=1)
    end
  end

  # Plot the sensitivity functions of interest
  pSigma = plot(f, valS, xscale = :log10, yscale = :log10, color = :red, w=1, label="\$\\sigma(S(j\\omega))\$")
  if size(valS,2) > 1
    for ii = 2:size(valS,2)
      plot(f, valS[:,ii], xscale = :log10, yscale = :log10, color = :red, w=1)
    end
  end
  plot!(f, valKS, xscale = :log10, yscale = :log10, color = :green, w=1, label="\$\\sigma(C(j\\omega)S(j\\omega)))\$")
  if size(valKS,2) > 1
    for ii = 2:size(valKS,2)
      plot(f, valKS[:,ii], xscale = :log10, yscale = :log10, color = :green, w=1)
    end
  end
  plot!(f, valT, xscale = :log10, yscale = :log10, color = :blue, w=1, label="\$\\sigma(T(j\\omega)))\$")
  if size(valT,2) > 1
    for ii = 2:size(valT,2)
      plot(f, valT[:,ii], xscale = :log10, yscale = :log10, color = :blue, w=1)
    end
  end

  # Visualize the weighting functions
  if isa(WS, LTISystem) || isa(WS, Number)
    valiWS = sigma(gamma/WS[1,1], f)[1]
    plot!(f, valiWS, xscale = :log10, yscale = :log10, color = :red, w=2, style=:dot, label="\$\\gamma\\sigma(W_S(j\\omega)^{-1})\$")
  end
  if isa(WU, LTISystem) || isa(WU, Number)
    valiWU = sigma(gamma/WU[1,1], f)[1]
    plot!(f, valiWU, xscale = :log10, yscale = :log10, color = :green, w=2, style=:dot, label="\$\\gamma\\sigma(W_U(j\\omega)^{-1})\$")
  end
  if isa(WT, LTISystem) || isa(WT, Number)
    valiWT = sigma(gamma/WT[1,1], f)[1]
    plot!(f, valiWT, xscale = :log10, yscale = :log10, color = :blue, w=2, style=:dot, label="\$\\gamma\\sigma(W_T(j\\omega)^{-1})\$")
  end

  # Visualize the step response
  time = [ii for ii in range(0,stop=tmax,length=1001)]
  stepy, stept, stepx = step(T, time)
  pStep1 = plot(stept, stepy, color = :blue, w=2)

  #stepy, stept, stepx = step(KS, tmax)
  #pStep2 = plot(stept[2:end], stepy[2:end], color = :red, w=2)

  l = @layout [ a ; b ; c ]
  plt=plot(pGain, pSigma, pStep1, layout=l, size=(600,800))
  gui()
  if isa(filename, String); savefig(filename); end
end
