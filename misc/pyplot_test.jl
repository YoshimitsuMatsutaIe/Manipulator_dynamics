# using PyPlot
# pygui(true)


using PyCall
pygui(true)
using PyPlot
# using PyCall
# plt = pyimport(matplotlib.pyplot)

# using PyCall
# @pyimport matplotlib.pyplot as plt

x = range(0, 2π, step=0.1)
y = []

# yの値をsin(x)とする
for i in x
    push!(y, sin(i))
end

# グラフ描画(x-y座標上にプロット)
#fig = plot(x, y)
plot(x,y)
show()
#savefig("temp.png")