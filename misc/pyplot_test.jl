#using PyPlot
using PyCall
@pyimport matplotlib.pyplot as plt


x = range(0, 2π, step=0.1)
y = []

# yの値をsin(x)とする
for i in x
    push!(y, sin(i))
end

# グラフ描画(x-y座標上にプロット)
#fig = plot(x, y)
fig = plt.plot(x,y)
plt.show()