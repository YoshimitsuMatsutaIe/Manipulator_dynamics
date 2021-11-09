using Plots

using ForwardDiff
eta_d = 1.0
eta_e = 0.5

pot(eta, e) = 1/eta * log(exp(eta*e) + exp(-eta*e))


e(y, y0) = abs(y0 - y)





yd = 0.5
ye = 1.0

potd(y) = pot(eta_d, e(y, yd))
pote(y) = pot(eta_e, e(y, ye))

pot2(y) = (pot(eta_d, e(y, yd)) + pot(eta_e, e(y, ye)))/2


# plot(potd, label="d")
# plot!(pote, label="e")
# plot!(pot2, label="sum")


nabla_pot2(y) = ForwardDiff.derivative(pot2, y)


y = -2.0:0.01:2.0
p2 = pot2.(y)
pd = potd.(y)
pe = pote.(y)

dp2 = [nabla_pot2(_y) for _y in y]

plot(y, p2, label="p2")
plot!(y, dp2,label="dp2")
plot!(y, pd, label="pd")
plot!(y, pe, label="pe")
