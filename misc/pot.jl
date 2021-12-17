using Plots

using ForwardDiff

x1 = 0  # 外力平衡位置
x2 = 0.5  # 目標位置（物体表面）



# a = 20
# b = x2^2 / 4
# c = 0


# function df(x)
#     return a*((x - x2/2)^2-b)^2
# end

# function f(x)
#     return 1/5*a*(x-x2/2)^5 -2/3*a*b*(x-x2/2)^3 + (a*b^2)*(x-x2/2) - offset(x)
# end

# function offset(x)
#     1/5*a*(-x2/2)^5 -2/3*a*b*(-x2/2)^3 + (a*b^2)*(-x2/2)
# end

# x = -0.5:0.01:5

# fs = f.(x)
# dfs = df.(x)


# plot(x, dfs, label="grad_pot", aspect_ratio=1, xlim=(-0.5,4), ylim=(-1,2))
# plot!(x, fs, label="pot")



# x1 = 0.01
# x2 = 0.1
# x3 = 4
# p2 = 0.01

# a = 3*p2/(5*x1*x2^4 - 50*x1*x2^3*x3 - 2*x2^5 + 5*x2^4*x3)
# b = (-5*p2*x1 - 5*p2*x2 - 5*p2*x3)/(5*x1*x2^4 - 50*x1*x2^3*x3 - 2*x2^5 + 5*x2^4*x3)
# c = (10*p2*x1*x2 + 10*p2*x1*x3 + 10*p2*x2*x3)/(5*x1*x2^4 - 50*x1*x2^3*x3 - 2*x2^5 + 5*x2^4*x3)
# d = -60*p2*x1*x3/(5*x1*x2^3 - 50*x1*x2^2*x3 - 2*x2^4 + 5*x2^3*x3)
# e = 0
# f = 0

# x = -1:0.01:4
# p(x) = a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f

# plot(x, p.(x), xlim=(-1, 4), ylim=(-1,1),label="pot")#, aspect_ratio=1)


xe = 0
xd = 0.3

eta=1
a=10

sig(x) = 1/(1+exp(-a*x))

dif_sig(x) = (1-sig(x))*sig(x)

pot(x) = 1/eta * log(exp(eta*x) + exp(-eta*x))

pot_2(x) = sig(x-xd)*pot(x-xd) + (1-sig(x-xd))*pot(x-xe)

plot(pot_2, label="pot", xlim=(0.0,3))
plot!(x -> ForwardDiff.derivative(pot_2, x), label="grad")
plot!(x -> sig(x-xd), label="sigmoid")
plot!(x -> dif_sig(x-xd), label="dif_sig")

