module TestPlots
using CustomTest
using ControlSystems, Plots
using VisualRegressionTests, ControlExamplePlots

default(show=false)

funcs, refs = getexamples()
res = genplots(funcs, refs, popup=false)

##Explicit enumeration for simpler debugging

#"bode.png"
println(res[1])
println(refs[1])
println(funcs[1])
@test  res[1] |> success
#"nyquist.png"
println(res[2])
println(refs[2])
println(funcs[2])
@test  res[2] |> success
#"sigma.png"
println(res[3])
println(refs[3])
println(funcs[3])
@test  res[3] |> success
#"nichols.png"
println(res[4])
println(refs[4])
println(funcs[4])
@test  res[4] |> success
#"step.png"
println(res[5])
println(refs[5])
println(funcs[5])
@test  res[5] |> success
#"impulse.png"
println(res[6])
println(refs[6])
println(funcs[6])
@test  res[6] |> success
#"lsim.png"
println(res[7])
println(refs[7])
println(funcs[7])
@test  res[7] |> success
#"margin.png"
println(res[8])
println(refs[8])
println(funcs[8])
@test  res[8] |> success
#"gangoffour.png"
println(res[9])
println(refs[9])
println(funcs[9])
@test  res[9] |> success
#"pzmap.png"
println(res[10])
println(refs[10])
println(funcs[10])
@test  res[10] |> success

end
