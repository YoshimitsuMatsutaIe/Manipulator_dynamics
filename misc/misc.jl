# using .hoge

# hoge.hello()

using Symbolics
@variables u[1:3]
function f(u)
  [u[1]-u[3],u[1]^2-u[2],u[3]+u[2]]
end
ex1, ex2 = build_function(f(u),u)
#write("function.jl", string(ex2))

using Latexify
latexify(ex2)