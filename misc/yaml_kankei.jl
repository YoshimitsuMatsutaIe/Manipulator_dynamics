using YAML
using Parameters





params = YAML.load_file("./config/sice.yaml")


@with_kw struct Env
    isLeft
    isWithMass
    TIME_SPAN
    TIME_INTERVAL
end


keytosymbol(x) = Dict(Symbol(k) => v for (k, v) in pairs(x))

e = params["sim_param"]
e_ = keytosymbol(e)


env = Env(;e_...)