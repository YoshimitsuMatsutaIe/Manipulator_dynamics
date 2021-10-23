# using PyCall
# pygui()
# popdisplay()
# using PyPlot

# #anim = pyimport("matplotlib.animation")
# @pyimport matplotlib.animation as anim

# x = range(0, 2π, step=0.1)
# y = sin.(x)
# z = cos.(x)



# function update(i)
#     println(i)
#     ax.cla()
#     ax.plot(x[i:i+4], y[i:i+4], z[i:i+4])
# end

# fig = plt.figure()
# ax = fig.add_subplot(projection="3d")
# ax.plot(x, y, z)

# ani = anim.FuncAnimation(fig, update, blit=false, interval=10, frames=10, repeat=false)

# show()




# using PyPlot
# using PyCall
# pygui()
# @pyimport matplotlib.animation as anim

# #Construct Figure and Plot Data
# fig = figure("MyFigure",figsize=(5,5))
# ax = fig.add_subplot(xlim = (0,10),ylim=(0,10))
# global line1 = ax[:plot]([],[],"r-")[1]
# global line2 = ax[:plot]([],[],"g-")[1]
# global line3 = ax[:plot]([],[],"b-")[1]

# # Define the init function, which draws the first frame (empty, in this case)
# function init()
#     global line1
#     global line2
#     global line3
#     line1[:set_data]([],[])
#     line2[:set_data]([],[])
#     line3[:set_data]([],[])
#     return (line1,line2,line3,Union{})  # Union{} is the new word for None
# end

# # Animate draws the i-th frame, where i starts at i=0 as in Python.
# function animate(i)
#     global line1
#     global line2
#     global line3
#     x = (0:i)/10.0
#     line1[:set_data](x,x)
#     line2[:set_data](1+x,x)
#     line3[:set_data](2+x,x)
#     return (line1,line2,line3,Union{})
# end

# # Create the animation object by calling the Python function FuncAnimaton
# myanim = anim.FuncAnimation(fig, animate, init_func=init, frames=100, interval=20)


# show()

using PyCall
popdisplay()
@pyimport matplotlib.animation as anim
using PyPlot


fig = figure(figsize=(4,4))
#ax = axes()

x = [0:0.01:2pi;]

# i=0,1,...,frames-1
function animate(i)
    println(i)
    clf()
    plot(sin.(x+i/10.0))
end

function init()
    plot(x,sin.(x))
end


myanim = anim.FuncAnimation(fig, animate, frames=200, init_func=init, interval=25, blit=false)
#myanim[:save]("test3.gif", writer="pillow")

show()
#showanim("test3.mp4")