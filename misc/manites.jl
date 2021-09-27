using Plots
A = [ x^2 + y^2 for x in -1:0.1:1, y in -1:0.2:1 ]  

@manipulate for θ in 0:2:90, h in 0:2:90
  wireframe(A, camera=(θ,h))
end
