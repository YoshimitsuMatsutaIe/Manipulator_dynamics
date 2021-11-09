using Plots


eta_d = 1.0
eta_e = 1.0

pot(eta, e) = 1/eta * log(exp(eta*e) + exp(-eta*e))


e(y, y0) = abs(y0 - y)





yd = 0.0
ye = 1.0

potd(y) = pot(eta_d, e(y, yd))
pote(y) = pot(eta_e, e(y, ye))

pot2(y) = (pot(eta_d, e(y, yd)) + pot(eta_e, e(y, ye)))/2


plot(potd, label="d")
plot!(pote, label="e")
plot!(pot2, label="sum")
