# 点でのrmpシミュレーション

using Dates
using LinearAlgebra
using Parameters
using Plots

include("rmp.jl")
include("utils.jl")




"""ルンゲクッタ法（4次）"""
function solve_RungeKutta(;f, x₀, t_span, Δt)

    t = range(t_span..., step = Δt)  # 時間軸
    x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列
    fg_norm = Vector{typeof(x₀[1])}(undef, length(t))
    fo_norm = Vector{typeof(x₀[1])}(undef, length(t))

    x[1] = x₀  # 初期値
    _, f_norms = ddx(x₀[1:2], x₀[3:4])
    fg_norm[1] = f_norms[1]
    fo_norm[1] = f_norms[2]

    for i in 1:length(x)-1
        
        k₁, _ = f(x[i])
        k₂, _ = f(x[i] .+ k₁ .* Δt/2)
        k₃, _ = f(x[i] .+ k₂ .* Δt/2)
        k₄, _ = f(x[i] .+ k₃ .* Δt)

        x[i+1] = x[i] .+ (k₁ .+ 2k₂ .+ 2k₃ .+k₄) .* Δt/6

        _, f_norms = ddx(x[i+1][1:2], x[i+1][3:4])
        fg_norm[i+1] = f_norms[1]
        fo_norm[i+1] = f_norms[2]
    end

    t, x, fg_norm, fo_norm
end







const p = (
    obs_avoidance = OriginalRMPCollisionAvoidance(
    scale_rep = 0.1,
    scale_damp = 0.1,
    ratio = 1.0,
    gain = 0.0,
    r = 0.1
    ),
    # obs_avoidance = RMPfromGDSCollisionAvoidance(
    #     rw = 0.5,
    #     sigma = 1.0,
    #     alpha = 1.0
    # ),
    # attractor = OriginalRMPAttractor(
    #     max_speed = 1.0,
    #     gain = 1.5,
    #     ddq_damp_r = 1.0,
    #     sigma_W = 1.0,
    #     sigma_H = 1.0,
    #     metric_damp_r = 0.1
    # )
    attractor = RMPfromGDSAttractor(
        max_speed = 1.0,
        gain = 1.5,
        f_alpha = 0.15,
        sigma_alpha = 1.0,
        sigma_gamma = 1.0,
        wu = 10.0,
        wl = 1.0,
        alpha = 0.15,
        epsilon = 1.0e-5
    )
)

const o = [0.0, 0.0]
const g = [-3.0, 3.0]


"""所望の加速度"""
function ddx(x, dx)

    fg, Mg = get_natural(p.attractor, x, dx, g)
    fo, Mo = get_natural(p.obs_avoidance, x, dx, o)

    #fo = zero(fg)
    #Mo = zero(Mg)

    root_f = fg .+ fo
    root_M = Mg .+ Mo

    z = pinv(root_M) * root_f
    return z, (norm(fg), norm(fo))
end


function dX(X)

    x = X[1:2]
    dx = X[3:4]
    _ddx, f_norms = ddx(x, dx)
    return [dx; _ddx], f_norms
end


function run()


    x0 = [2.5, -3.2]
    #dx0 = [-1.0, 1.0]

    #x0 = [-3.0, 3.0]
    dx0 = [0.0, 0.0]

    t, X, fg, fo = solve_RungeKutta(
        f = dX,
        x₀ = [x0; dx0],
        t_span = (0.0, 400.0),
        Δt = 0.01
    )


    # グラフ化
    x, y, dx, dy = split_vec_of_arrays(X)
    
    fig = plot(x, y, label = "x", legend=:outerright)
    scatter!(fig, [x[1]], [y[1]], label="x0")
    scatter!(fig, [x[end]], [y[end]], label="xend")
    scatter!(fig, [g[1]], [g[2]], markershape=:star6, label="g")
    scatter!(fig, [o[1]], [o[2]], label="o")

    fig2 = plot(t, x, label="_x", legend=:outerright)
    plot!(fig2, t, y, label="_y")
    plot!(fig2, t, dx, label="dx")
    plot!(fig2, t, dy, label="dy")

    fig3 = plot(t, fg, label="fg", legend=:outerright)
    plot!(fig3, t, fo, label="fo")

    plot(
        fig, fig2, fig3,
        layout=(3, 1), size=(500,900)
    )
end






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


path = get_time_string()  # 実行時のデータ保存パス

@time run()

