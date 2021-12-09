"""siceのシミュレーション"""

using Plots
using StaticArrays
using ArraysOfArrays
using Parameters
using LinearAlgebra
using Random

include("./sice_kinematics.jl")
include("./sice_robot_dynamics.jl")
include("./rmp.jl")

using .SiceKinematics
using .SiceDynamics


const q_min = [(-5/4)π+(1/2)π, (-5/4)π, (-5/4)π, (-5/4)π]
const q_max = [(5/4)π+(1/2)π, (5/4)π, (5/4)π, (5/4)π]
const q_neutral = [(1/2)π, 0.0, 0.0, 0.0]


function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end

function calc_min_dis_to_obs(xs::Vector{Vector{T}}, xos::Vector{Vector{T}}) where T

    d = Matrix{T}(undef, length(xs), length(xos))
    for (i, x) in enumerate(xs)
        for (j, xo) in enumerate(xos)
            d[i, j] = norm(x - xo)
        end
    end
    
    return minimum(d)
end

"""ジョイント制限を守れてるかチェック"""
function check_JointLimitation(q::Vector{T}) where T
    return (q .< q_min) .| (q_max .< q)
end

"""円物体が与える外力"""
function Fc_from_cicle(
    x::Vector{T}, dx::Vector{T},
    center::Vector{T}, r::T,
    K::T, D::T
    ) where T
    s = norm(x .- center)

    if s >= r
        return zero(x)
    else
        return @. K * (x - center) + D * dx
    end
end




"""実験データ"""
struct Data{T}
    t::StepRangeLen{T}
    q::Vector{Vector{T}}
    dq::Vector{Vector{T}}
    ddq::Vector{Vector{T}}
    desired_ddq::Vector{Vector{T}}
    u::Vector{Vector{T}}
    x1::Vector{Vector{T}}
    x2::Vector{Vector{T}}
    x3::Vector{Vector{T}}
    x4::Vector{Vector{T}}
    dx1::Vector{Vector{T}}
    dx2::Vector{Vector{T}}
    dx3::Vector{Vector{T}}
    dx4::Vector{Vector{T}}
    J1::Vector{Matrix{T}}
    J2::Vector{Matrix{T}}
    J3::Vector{Matrix{T}}
    J4::Vector{Matrix{T}}
    dJ1::Vector{Matrix{T}}
    dJ2::Vector{Matrix{T}}
    dJ3::Vector{Matrix{T}}
    dJ4::Vector{Matrix{T}}
    f0::Vector{Vector{T}}
    f1::Vector{Vector{T}}
    f2::Vector{Vector{T}}
    f3::Vector{Vector{T}}
    f4::Vector{Vector{T}}
    M0::Vector{Matrix{T}}
    M1::Vector{Matrix{T}}
    M2::Vector{Matrix{T}}
    M3::Vector{Matrix{T}}
    M4::Vector{Matrix{T}}
    error::Vector{T}
    min_dit_to_obs::Vector{T}
    jl::Vector{BitVector}
    F_distur::Vector{Vector{T}}  # トルクに入る外乱
    Fc::Vector{Vector{T}}  # 対象に与えた力（i.e. 対象が受ける力 = -Fc）
end



"""全部実行"""
function run_simulation(TIME_SPAN::T=1000.0, Δt::T=0.01, isImpedance::Bool=true) where T

    # 初期値
    q₀ = q_neutral
    dq₀ = zero(q₀)
    ddq₀ = zero(q₀)
    desired_ddq₀ = zero(q₀)
    u₀ = zero(q₀)
    F_distur₀ = zero(q₀)

    # 目標値
    xd = [2.0, 1.0]
    #xd = [2.0, 1.0].- [0.0, 0.01]

    # コンプライアンス中心
    xe = xd .- [0.0, 0.005]

    # 障害物
    xo = [
        [1.9, 2.0],
        [2.1, 2.0]
    ] * 1000

    # 対象物
    # box_l = 2.0
    # box_h = 1.0
    # box_center = [2.0, 0.5]

    circle = (
        r = 0.5, x = 2.0, y = 0.5, K = 1.0, D = 10.0,
    )  # 円の物体の情報




    # rmpのパラメータ


    obs_avoidance = RMPfromGDSCollisionAvoidance(
        rw = 0.3,
        sigma = 1.0,
        alpha = 2.0,
    )

    jl_avoidance = RMPfromGDSJointLimitAvoidance(
        gamma_p = 0.05,
        gamma_d = 0.01,
        lambda = 1.0,
        sigma = 0.1,
    )

    if isImpedance
        # インピーダンス特性の決定
        zeta_d = 0.9  # 所望の減衰係数
        omega_d = 10.0
        dd = 1.0

        de = circle.D
        ke = circle.K

        md = (de + dd)/(2*omega_d*zeta_d)
        kd = (-ke*zeta_d + omega_d*(de + dd)/2)/zeta_d

        println("md = ", md)
        println("kd = ", kd)
        println("dd = ", dd)


        impedance = RMPfromGDSImpedance(
            M_d = Matrix{T}(I, 2, 2) * md,
            D_d = Matrix{T}(I, 2, 2) * kd,
            P_d = Matrix{T}(I, 2, 2) * kd,
            D_e = Matrix{T}(I, 2, 2) * kd,
            P_e = Matrix{T}(I, 2, 2) * ke,
            a=10.0,
            eta_d=1.0,
            eta_e=1.0,
            f_alpha=0.15,
            sigma_alpha=1.0,
            sigma_gamma=1.0,
            wu=5.0,
            wl=1.0,
            alpha=0.15,
            epsilon=0.05,
        )

    else
        attractor = RMPfromGDSAttractor(
            max_speed = 2.0,
            gain = 10.0,
            f_alpha = 0.15,
            sigma_alpha = 1.0,
            sigma_gamma = 1.0,
            wu = 5.0,
            wl = 1.0,
            alpha = 0.15,
            epsilon = 0.05,
        )

    end



    t = range(0.0, TIME_SPAN, step=Δt)  # 時間列

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{BitVector}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
    )



    # 初期値代入
    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = ddq₀
    data.desired_ddq[1] = desired_ddq₀
    data.u[1] = u₀

    data.x1[1], data.x2[1], data.x3[1], data.x4[1] = calc_x(q₀)
    data.J1[1], data.J2[1], data.J3[1], data.J4[1] = calc_J(q₀)
    data.dJ1[1], data.dJ2[1], data.dJ3[1], data.dJ4[1] = calc_dJ(q₀, dq₀)

    data.dx1[1] = data.J1[1] * data.dq[1]
    data.dx2[1] = data.J2[1] * data.dq[1]
    data.dx3[1] = data.J3[1] * data.dq[1]
    data.dx4[1] = data.J4[1] * data.dq[1]

    data.f0[1] = zero(q₀)
    data.f1[1] = zeros(T, 2)
    data.f2[1] = zeros(T, 2)
    data.f3[1] = zeros(T, 2)
    data.f4[1] = zeros(T, 2)

    data.M0[1] = zeros(T, 4, 4)
    data.M1[1] = zeros(T, 2, 2)
    data.M2[1] = zeros(T, 2, 2)
    data.M3[1] = zeros(T, 2, 2)
    data.M4[1] = zeros(T, 2, 2)

    data.error[1] = norm(data.x4[1] - xd)
    data.min_dit_to_obs[1] = calc_min_dis_to_obs(
        [data.x1[1], data.x2[1], data.x3[1], data.x4[1]],
        xo
    )
    data.jl[1] = check_JointLimitation(q₀)
    data.F_distur[1] = zero(u₀)
    data.Fc[1] = Fc_from_cicle(
        data.x4[1], data.dx4[1],
        [circle.x, circle.y], circle.r, circle.K, circle.D
    )


    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #println("t = ", data.t[i])


        # 状態を更新
        data.q[i+1] = data.q[i] .+ data.dq[i] .* Δt
        data.dq[i+1] = data.dq[i] .+ data.ddq[i] .*Δt

        data.x1[i+1], data.x2[i+1], data.x3[i+1], data.x4[i+1] = calc_x(data.q[i+1])
        data.J1[i+1], data.J2[i+1], data.J3[i+1], data.J4[i+1] = calc_J(data.q[i+1])
        data.dJ1[i+1], data.dJ2[i+1], data.dJ3[i+1], data.dJ4[i+1] = calc_dJ(data.q[i+1], data.dq[i+1])

        data.dx1[i+1] = data.J1[i+1] * data.dq[i+1]
        data.dx2[i+1] = data.J2[i+1] * data.dq[i+1]
        data.dx3[i+1] = data.J3[i+1] * data.dq[i+1]
        data.dx4[i+1] = data.J4[i+1] * data.dq[i+1]



        # rmpを計算

        # rmpをセット
        # root
        data.f0[i+1], data.M0[i+1] = get_natural(
            jl_avoidance, data.q[i+1], data.dq[i+1], q_max, q_min, q_neutral
        )
        data.f1[i+1] = zeros(T, 2)
        data.f2[i+1] = zeros(T, 2)
        data.f3[i+1] = zeros(T, 2)
        data.f4[i+1] = zeros(T, 2)
        data.M1[i+1] = zeros(T, 2, 2)
        data.M2[i+1] = zeros(T, 2, 2)
        data.M3[i+1] = zeros(T, 2, 2)
        data.M4[i+1] = zeros(T, 2, 2)


        for o in xo
            f, M = get_natural(obs_avoidance, data.x1[i+1], data.dx1[i+1], o)
            @. data.f1[i+1] += f
            @. data.M1[i+1] += M

            f, M = get_natural(obs_avoidance, data.x2[i+1], data.dx2[i+1], o)
            @. data.f2[i+1] += f
            @. data.M2[i+1] += M

            f, M = get_natural(obs_avoidance, data.x3[i+1], data.dx3[i+1], o)
            @. data.f3[i+1] += f
            @. data.M3[i+1] += M

            f, M = get_natural(obs_avoidance, data.x4[i+1], data.dx4[i+1], o)
            @. data.f4[i+1] += f
            @. data.M4[i+1] += M
        end

        if isImpedance
            f, M = get_natural(impedance, data.x4[i+1], xd, xe, data.dx4[i+1])  # インピーダンス
        else
            f, M = get_natural(attractor, data.x4[i+1], data.dx4[i+1], xd)  # アトラクタ
        end
        @. data.f4[i+1] += f
        @. data.M4[i+1] += M
        #println("attractor = ", f)


        # #println("attractor = ", f)

        # pullback演算
        root_f = data.f0[i+1]
        root_M = data.M0[i+1]

        _rf, _rM = pullbacked_rmp(data.f1[i+1], data.M1[i+1], data.J1[i+1], data.dJ1[i+1], data.dq[i+1])
        @. root_f += _rf
        @. root_M += _rM

        _rf, _rM = pullbacked_rmp(data.f2[i+1], data.M2[i+1], data.J2[i+1], data.dJ2[i+1], data.dq[i+1])
        @. root_f += _rf
        @. root_M += _rM

        _rf, _rM = pullbacked_rmp(data.f3[i+1], data.M3[i+1], data.J3[i+1], data.dJ3[i+1], data.dq[i+1])
        @. root_f += _rf
        @. root_M += _rM

        _rf, _rM = pullbacked_rmp(data.f4[i+1], data.M4[i+1], data.J4[i+1], data.dJ4[i+1], data.dq[i+1])
        @. root_f += _rf
        @. root_M += _rM

        # resolve演算
        data.desired_ddq[i+1] = pinv(root_M) * root_f

        data.u[i+1] = calc_torque(data.q[i+1], data.dq[i+1], data.desired_ddq[i+1])

        data.F_distur[i+1] = zero(u₀)  # 外乱無し
        #data.F_distur[i+1] = rand(T, 4)*0.1

        data.Fc[i+1] = Fc_from_cicle(
            data.x4[i+1], data.dx4[i+1],
            [circle.x, circle.y], circle.r, circle.K, circle.D
        )

        data.ddq[i+1] = calc_real_ddq(
            u = data.u[i+1],
            q = data.q[i+1],
            dq = data.dq[i+1],
            F = data.F_distur[i+1],
            Fc = data.Fc[i+1],
            Jend = data.J4[i+1],
        )

        data.error[i+1] = norm(data.x4[i+1] - xd)
        #println("error = ", data.error[i+1])
        data.min_dit_to_obs[i+1] = calc_min_dis_to_obs(
            [data.x1[i+1], data.x2[i+1], data.x3[i+1], data.x4[i+1]],
            xo,
        )
        data.jl[i+1] = check_JointLimitation(data.q[i+1])




    end



    # グラフ化

    x, y, z, w = split_vec_of_arrays(data.q)
    fig_q = plot(xlim=(0, TIME_SPAN), legend=:outerright, ylabel="q")
    plot!(fig_q, data.t, x, label="_q1")
    plot!(fig_q, data.t, y, label="_q2")
    plot!(fig_q, data.t, z, label="_q3")
    plot!(fig_q, data.t, w, label="_q4")

    x, y, z, w = split_vec_of_arrays(data.dq)
    fig_dq = plot(xlim=(0, TIME_SPAN), legend=:outerright, ylabel="dq")
    plot!(fig_dq, data.t, x, label="_w1")
    plot!(fig_dq, data.t, y, label="_w2")
    plot!(fig_dq, data.t, z, label="_w3")
    plot!(fig_dq, data.t, w, label="_w4")

    x, y, z, w = split_vec_of_arrays(data.ddq)
    fig_ddq = plot(xlim=(0, TIME_SPAN), legend=:outerright, ylabel="ddq")
    plot!(fig_ddq, data.t, x, label="_a1")
    plot!(fig_ddq, data.t, y, label="_a2")
    plot!(fig_ddq, data.t, z, label="_a3")
    plot!(fig_ddq, data.t, w, label="_a4")

    x, y, z, w = split_vec_of_arrays(data.desired_ddq)
    fig_desired_ddq = plot(xlim=(0, TIME_SPAN), legend=:outerright, ylabel="desired_ddq")
    plot!(fig_desired_ddq, data.t, x, label="a*1")
    plot!(fig_desired_ddq, data.t, y, label="a*2")
    plot!(fig_desired_ddq, data.t, z, label="a*3")
    plot!(fig_desired_ddq, data.t, w, label="a*4")

    x, y, z, w = split_vec_of_arrays(data.u)
    fig_u = plot(xlim=(0, TIME_SPAN), legend=:outerright, ylabel="u")
    plot!(fig_u, data.t, x, label="_u1")
    plot!(fig_u, data.t, y, label="_u2")
    plot!(fig_u, data.t, z, label="_u3")
    plot!(fig_u, data.t, w, label="_u4")


    f0 = norm.(data.f0)
    f1 = norm.(data.f1)
    f2 = norm.(data.f2)
    f3 = norm.(data.f3)
    f4 = norm.(data.f4)
    _fs = [0.0]
    f_max = maximum(append!(_fs, f0, f1, f2, f3, f4))
    fig_f = plot(
        xlim=(0, TIME_SPAN), ylim=(0, f_max),
        legend=:outerright, ylabel="f"
    )
    plot!(fig_f, data.t, f0, label="roo")
    plot!(fig_f, data.t, f1, label="_f1")
    plot!(fig_f, data.t, f2, label="_f2")
    plot!(fig_f, data.t, f3, label="_f3")
    plot!(fig_f, data.t, f4, label="_ee")

    fig_error = plot(
        data.t, data.error,
        label="err", ylabel="error",
        xlim=(0, TIME_SPAN), ylim=(0, maximum(data.error)),
        legend=:outerright,
    )

    # # 物体に与えた力をノルムでplot
    # Fcs = norm.(data.Fc)
    # fig_Fc = plot(
    #     data.t, Fcs,
    #     label="_Fc", ylabel="Fc",
    #     xlim=(0, TIME_SPAN), ylim=(0, maximum(Fcs)),
    #     legend=:outerright,
    # )

    println(data.Fc[3])
    x, y = split_vec_of_arrays(data.Fc)
    fig_Fc = plot(
        label="_Fc", ylabel="Fc",
        xlim=(0, TIME_SPAN),
        legend=:outerright,
    )
    plot!(fig_Fc, data.t, x, label="Fcx")
    plot!(fig_Fc, data.t, y, label="Fcy")

    fig_dis_to_obs = plot(
        data.t, data.min_dit_to_obs,
        label="obs", ylabel="min dis to obs",
        xlim=(0, TIME_SPAN), ylim=(0, maximum(data.min_dit_to_obs)),
        legend=:outerright,
    )

    fig = plot(
        fig_q, fig_dq, fig_ddq,
        fig_desired_ddq, fig_u, fig_f,
        fig_error, fig_Fc, fig_dis_to_obs,
        layout=(9, 1),
        size=(500, 2200)
    )

    if isImpedance
        name = "sice_simple_proposed.png"
    else
        name = "sice_simple_conventional.png"
    end

    savefig(fig, name)


    # アニメ作成

    # 少し準備
    """箱の形"""
    function rectangle(w, h, x, y)
        Shape(x .+ [-w/2,w/2,w/2,-w/2], y .+ [-h/2,-h/2,h/2,h/2])
    end

    """楕円"""
    function ellipse(x, y, w, h)
        ang = range(0, 2π, length = 60)
        Shape(w*sin.(ang).+x, h*cos.(ang).+y)
    end


    """1フレームを描写"""
    function draw_frame(i)
        arm = [
            [0.0, 0.0],
            data.x1[i],
            data.x2[i],
            data.x3[i],
            data.x4[i]
        ]

        x, y = split_vec_of_arrays(arm)
        fig = plot(
            x, y,
            marker=:circle,
            aspect_ratio = 1,
            xlabel = "X[m]", ylabel = "Y[m]",
        )

        scatter!(
            fig,
            [xd[1]], [xd[2]],
            markershape=:star6,
        )

        x, y = split_vec_of_arrays(xo)
        scatter!(
            fig,
            x, y
        )

        # 対象物
        #plot!(rectangle(box_l, box_h, box_center[1], box_center[2]), opacity=.5)
        plot!(ellipse(circle.x, circle.y, circle.r, circle.r), opacity=.5)

        x_max = 4.1
        x_min = -4.1
        y_max = 4.1
        y_min = -4.1
        max_range = max(x_max-x_min, y_max-y_min)*0.5
        x_mid = (x_max + x_min) / 2
        y_mid = (y_max + y_min) / 2

    
        plot!(
            fig,
            xlims=(x_mid-max_range, x_mid+max_range),
            ylims=(y_mid-max_range, y_mid+max_range),
            legend = false,
            size=(600, 600),
        )


        # rmpも描画

        scale = 10.0

        f = [data.x1[i], data.x1[i] .+ data.f1[i] .* scale]
        x, y = split_vec_of_arrays(f)
        plot!(
            fig,
            x, y
        )

        f = [data.x2[i], data.x2[i] .+ data.f2[i] .* scale]
        x, y = split_vec_of_arrays(f)
        plot!(
            fig,
            x, y
        )

        f = [data.x3[i], data.x3[i] .+ data.f3[i] .* scale]
        x, y = split_vec_of_arrays(f)
        plot!(
            fig,
            x, y
        )

        f = [data.x4[i], data.x4[i] .+ data.f4[i] .* scale]
        x, y = split_vec_of_arrays(f)
        plot!(
            fig,
            x, y
        )

        # function make_metric_cicle(i, M, q_global)
        #     e = eigen(M)
        #     axis_len = 1 / sqrt.(e.Vector) * 0.1
        #     max_len = maximum(axis_len)
        #     scale = minimum(2/max_len, 1)
        #     M_xis_len = axis_len .* scale
        #     angle = atan(e.Vector[2], e.Vector[1]) + q_global |>
        #     rad2deg

        # end

        id = findall(x->x==true, data.jl[i])
        plot!(fig, title = string(round(data.t[i], digits=2)) * "[s]" * ", " * string(id))
        return fig
    end

    println("アニメ作成中...")
    # 枚数決める
    #println(data.t)
    epoch_max = 100
    epoch = length(data.t)
    if epoch < epoch_max
        step = 1
    else
        step = div(epoch, epoch_max)
    end

    #println(step)

    anim = Animation()
    @gif for i in 1:step:length(data.q)
        _fig = draw_frame(i)
        frame(anim, _fig)
    end


    if isImpedance
        name = "sice_animation_proposed.gif"
    else
        name = "sice_animation_conventional.gif"
    end

    gif(anim, name, fps = 60)
    
    println("アニメ作成完了")

    return data
end



println("実行中...")
@time data = run_simulation()
println("実行終了!!")