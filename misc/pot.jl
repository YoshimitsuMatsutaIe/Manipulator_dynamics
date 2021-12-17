using Plots



x1 = 0  # 外力平衡位置
x2 = 0.5  # 目標位置（物体表面）



a = 20
b = x2^2 / 4
c = 0


function df(x)
    return a*((x - x2/2)^2-b)^2
end

function f(x)
    return 1/5*a*(x-x2/2)^5 -2/3*a*b*(x-x2/2)^3 + (a*b^2)*(x-x2/2) - offset(x)
end

function offset(x)
    1/5*a*(-x2/2)^5 -2/3*a*b*(-x2/2)^3 + (a*b^2)*(-x2/2)
end

x = -0.5:0.01:5

fs = f.(x)
dfs = df.(x)


plot(x, dfs, label="grad_pot", aspect_ratio=1, xlim=(-0.5,4), ylim=(-1,2))
plot!(x, fs, label="pot")
