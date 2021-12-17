using Plots

eta = 1
a = 10

sig(x) = 1/(1+exp(-a*x))
pot(x) = 1/eta * log(exp(eta*x)+exp(-eta*x))

xd = 2

pot_2(x) = (sig(x-xd))*pot(x-xd) + (1-sig(x-xd))*pot(x)


pot_3(x) = 



plot(xlim=(0,6), ylim=(0,5))
plot!(x -> pot(x), label="Pe")
plot!(x -> pot(x-xd), label="Pd")
plot!(x -> pot_2(x))
plot!(x -> sig(x-xd), label="sig")