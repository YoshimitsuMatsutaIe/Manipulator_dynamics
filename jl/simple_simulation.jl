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



"""multi-robot-rmpのrmp"""
function multi_attractor_rmp(x::Vector{T}, dx::Vector{T}, g::Vector{T}) where T
    wu = 10.0
    wl = 1.0
    sigma = 1.0
    alpha = 1.0
    eta = 2.0
    gain = 1.0
    tol = 0.005

    z = x - g
    dz = dx
    #J = Matrix{T}(I, 2, 2)
    #dJ = zeros(T, 2, 2)

    z_norm = norm(z)
    beta = exp(-z_norm^2 / 2 / (sigma^2))
    w = (wu - wl) * beta + wl
    s = (1-exp(-2*alpha*z_norm)) / (1+exp(-2*alpha*z_norm))

    G = Matrix{T}(I, 2, 2) * w
    
    if z_norm > tol
        grad_phi = s / z_norm * w * z * gain
    else
        grad_phi = 0.0
    end

    Bx_dot = eta * w * dz
    grad_w = -beta * (wu-wl) / sigma^2*z

    dz_norm = norm(dz)

    xi = -0.5 * (dz_norm^2 * grad_w - 2 * dz * dz' * dz * dz' * grad_w)

    M = G
    f = -grad_phi - Bx_dot - xi

    return f, M
end


"""multi-robot-rmpのrmp"""
function multi_avoidance_rmp(x::Vector{T}, dx::Vector{T}, o::Vector{T}) where T
    R = 1.0
    alpha = 1e-5
    eta = 0.0
    epsilon = 0.2

    z = norm(x-o)/R -1
    J = 1/norm(x-o) * (x-o)' / R
    dJ = (dx' * (-1/norm(x-o)^3 * (x-o) * (x-o)' + 1/norm(x-o) * Matrix{T}(I, 2, 2))) / R
    dz = J * dx


    if z < 0
        w = 1.0e10
        grad_w = 0.0
    else
        w = 1/z^4
        grad_w = -4/z^5
    end
    
    u = epsilon + min(0.0, dz) * dz
    g = w * u

    grad_u = 2 * min(0.0, dz)
    grad_phi = alpha * w * grad_w
    xi = 0.5 * dz^2 * u * grad_w

    _M = g + 0.5 * dz * w * grad_u
    M = min(max(_M, -1.0e5), 1.0e5)

    Bx_dot = eta * g * dz

    _f = grad_phi - xi -Bx_dot
    f = min(max(_f, -1.0e10), 1.0e10)

    pulled_f = J' * (f .- M * dJ * dx)
    pulled_M = J' * M * J

    return pulled_f, pulled_M
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

    #fg, Mg = get_natural(p.attractor, x, dx, g)
    #fo, Mo = get_natural(p.obs_avoidance, x, dx, o)

    #fo = zero(fg)
    #Mo = zero(Mg)

    fg, Mg = multi_attractor_rmp(x, dx, g)
    fo, Mo = multi_avoidance_rmp(x, dx, o)


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
    dx0 = [-1.0, 1.0]

    #x0 = [-3.0, 3.0]
    #dx0 = [0.0, 0.0]

    t, X, fg, fo = solve_RungeKutta(
        f = dX,
        x₀ = [x0; dx0],
        t_span = (0.0, 5.0),
        Δt = 0.001
    )


    # グラフ化
    x, y, dx, dy = split_vec_of_arrays(X)
    
    fig = plot(x, y, label = "x", legend=:outerright)
    scatter!(fig, [x[1]], [y[1]], label="x0")
    scatter!(fig, [x[end]], [y[end]], label="xend")
    scatter!(fig, [g[1]], [g[2]], markershape=:star6, label="g")
    scatter!(fig, [o[1]], [o[2]], label="o")

    function circle()
        theta = LinRange(0, 2*pi, 500)
        o[1] .+ 1*sin.(theta), o[2] .+ 1*cos.(theta)
    end
    plot!(
        fig, circle(), seriestype=[:shape,], lw=0.5,
        legend=false, fillalpha=0.2, aspect_ratio=1,
    
    )

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

