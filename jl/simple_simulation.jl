# 点でのrmpシミュレーション

using Dates
using LinearAlgebra
using Parameters
using Plots

include("rmp.jl")
include("utils.jl")


"""データ保存先のパス

シミュレーション結果を保存するディレクトリのパスを返す  
ないなら新規作成  
"""
function get_time_string()
    _nd = now()  # 時刻取得
    _Y = _nd |> Dates.year |> string
    _M = _nd |> Dates.month |> string
    _D = _nd |> Dates.day |> string
    _h = _nd |> Dates.hour |> string
    _m = _nd |> Dates.minute |> string
    _s = _nd |> Dates.second |> string
    t = _Y * _M * _D * "-" * _h * _m * _s

    linpath = "../result_of_manipulator_dynamics/" * _Y * _M * _D * "/"
    if !ispath(linpath)
        mkdir(linpath)
    end
    path = linpath * t
    return path
end


function solve_RungeKutta(;f, x₀, t_span, Δt)
    """ルンゲクッタ法（4次）"""

    t = range(t_span..., step = Δt)  # 時間軸
    x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

    x[1] = x₀  # 初期値
    for i in 1:length(x)-1
        
        k₁ = f(t[i], x[i])
        k₂ = f(t[i]+Δt/2, x[i] .+ k₁ .* Δt/2)
        k₃ = f(t[i]+Δt/2, x[i] .+k₂ .* Δt/2)
        k₄ = f(t[i]+Δt, x[i] .+ k₃ .* Δt)
        
        x[i+1] = x[i] .+ (k₁ .+ 2k₂ .+ 2k₃ .+k₄) .* Δt/6
    end

    t, x
end

const p = (
    obs_avoidance = RMPfromGDSCollisionAvoidance(
        rw = 0.5,
        sigma = 1.0,
        alpha = 0.1
    ),
    attractor = RMPfromGDSAttractor(
        max_speed = 1.0,
        gain = 1.0,
        f_alpha = 0.15,
        sigma_alpha = 2.0,
        sigma_gamma = 2.0,
        wu = 1.0,
        wl = 0.01,
        alpha = 0.15,
        epsilon = 1.0e-5
    )
)

const o = [1.2, 1.0]
const g = [2.0, 2.0]



"""所望の加速度"""
function ddx(x, dx, g, o, p)
    root_f = zero(x)
    root_M = zeros(typeof(x[1]), length(x), length(x))

    f, M = get_natural(p.attractor, x, dx, g)
    
    #println(f)
    #println(M)
    @. root_f += f
    @. root_M += M


    f, M = get_natural(p.obs_avoidance, x, dx, o)
    #println(f)
    #println(M)
    @. root_f += f
    @. root_M += M

    z = pinv(root_M) * root_f
    return z
end


function dX(t, X)

    x = X[1:2]
    dx = X[3:4]
    _ddx = ddx(x, dx, g, o, p)
    return [dx; _ddx]
end



function run()


    x0 = [0.0, 0.0]
    dx0 = [0.0, 0.5]


    t, X = solve_RungeKutta(
        f = dX,
        x₀ = [x0; dx0],
        t_span = (0.0, 100.0),
        Δt = 0.1
    )


    # グラフ化
    x, y, dx, dy = split_vec_of_arrays(X)
    
    fig = plot(x, y, label = "x", legend=:outerright)
    scatter!(fig, [g[1]], [g[2]], markershape=:star6, label="g")
    scatter!(fig, [o[1]], [o[2]], label="o")

end











path = get_time_string()  # 実行時のデータ保存パス

run()