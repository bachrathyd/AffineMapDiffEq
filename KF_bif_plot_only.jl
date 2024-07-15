##

using Plots
using LaTeXStrings
using MAT


file = matopen("q10_bif.mat")
xH = read(file, "xH")
yH = read(file, "yH")
xFup = read(file, "xFup")
yFup = read(file, "yFup")
xFdown = read(file, "xFdown")
yFdown = read(file, "yFdown")
xFix = read(file, "xFix")
yFix = read(file, "yFix")
xFup2 = read(file, "xFup2")
yFup2 = read(file, "yFup2")
xFdown2 = read(file, "xFdown2")
yFdown2 = read(file, "yFdown2")
# note that this does NOT introduce a variable ``varname`` into scope
close(file)





#plotly()
gr()
scatter([xH], [yH])
plot!(xFup[:], yFup[:])
plot!(xFdown[:], yFdown[:])
plot!(xFix[:], yFix[:])
plot!(xFup2[:], yFup2[:])
plot!(xFdown2[:], yFdown2[:])
plot!(xlabel=L"β", ylabel=L"y_{1}", legend=false)

savefig("bif_diag.svg")