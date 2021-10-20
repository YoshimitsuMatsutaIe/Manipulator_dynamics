using CPUTime
using Plots
using LinearAlgebra


include("../utils.jl")


"""DHパラメータ（修正？）"""
mutable struct DHparam{T}
    α::T
    a::T
    d::T
    θ::T
end


# パラメータ
const L = 278e-3
const h = 64e-3
const H = 1104e-3
const L0 = 270.35e-3
const L1 = 69e-3
const L2 = 364.35e-3
const L3 = 69e-3
const L4 = 374.29e-3
const L5 = 10e-3
const L6 = 368.3e-3

const q_neutral = [0.0; -31.0; 0.0; 43.0; 0.0; 72.0; 0.0] * pi/180  # ニュートラルの姿勢
const q_max = [51.0, 60.0, 173.0, 150.0, 175.0, 120.0, 175.0] * pi/180
const q_min = [-141.0, -123.0, -173.0, -3.0, -175.0, -90.0, -175.0] * pi/180

const DHparams_neutral = [
    DHparam(0.0, 0.0, 0.0, q_neutral[1])
    DHparam(-pi/2, L1, 0.0, q_neutral[2]+pi/2)
    DHparam(pi/2, 0.0, L2, q_neutral[3])
    DHparam(-pi/2, L3, 0.0, q_neutral[4])
    DHparam(pi/2, 0.0, L4, q_neutral[5])
    DHparam(-pi/2, L5, 0.0, q_neutral[6])
    DHparam(pi/2, 0.0, 0.0, q_neutral[7])
]

# 制御点
const cpoints_local = (
    (
        [0.0; L1/2; -L0/2; 1.0],
        [0.0; -L1/2; -L0/2; 1.0],
        [L1/2; 0.0; -L0/2; 1.0],
        [-L1/2; 0.0; -L0/2; 1.0],
    ),  # 1
    (
        [0.0; 0.0; L3/2; 1.0],
        [0.0; 0.0; -L3/2; 1.0],
    ),  # 2
    (
        [0.0; L3/2; -L2*2/3; 1.0],
        [0.0; -L3/2; -L2*2/3; 1.0],
        [L3/2; 0.0; -L2*2/3; 1.0],
        [-L3/2; 0.0; -L2*2/3; 1.0],
        [0.0; L3/2; -L2*1/3; 1.0],
        [0.0; -L3/2; -L2*1/3; 1.0],
        [L3/2; 0.0; -L2*1/3; 1.0],
        [-L3/2; 0.0; -L2*1/3; 1.0],
    ),  # 3
    (
        [0.0; 0.0; L3/2; 1.0],
        [0.0; 0.0; -L3/2; 1.0],
    ),  # 4
    (
        [0.0; L5/2; -L4/2; 1.0],
        [0.0; -L5/2; -L4/2; 1.0],
        [L5/2; 0.0; -L4/2; 1.0],
        [-L5/2; 0.0; -L4/2; 1.0],
    ),  # 5
    (
        [0.0; 0.0; L5/2; 1.0],
        [0.0; 0.0; -L5/2; 1.0]
    ),  # 6
    (
        [0.0; L5/2; L6/2; 1.0],
        [0.0; -L5/2; L6/2; 1.0],
        [L5/2; 0.0; L6/2; 1.0],
        [-L5/2; 0.0; L6/2; 1.0],
    ),  # 7
    (
        [0.0; 0.0; L6/3; 1.0],
    ),  # 8
)


# 右手
const HTM_BR_Wo = [
    -sqrt(2)/2 sqrt(2)/2 0.0 -L
    -sqrt(2)/2 -sqrt(2)/2 0.0 -h
    0.0 0.0 1.0 H
    0.0 0.0 0.0 1.0
]

const HTM_0_BR = [
    1.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 L0
    0.0 0.0 0.0 1.0
]

# 左手
const HTM_BL_Wo = [
    sqrt(2)/2 sqrt(2)/2 0.0 L
    -sqrt(2)/2 sqrt(2)/2 0.0 -h
    0.0 0.0 1.0 H
    0.0 0.0 0.0 1.0
]

const HTM_0_BL = HTM_0_BR


const HTM_GR_7 = [
    1.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 L6
    0.0 0.0 0.0 1.0
]

const HTM_A = [
    0.0 -1.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
]  # 偏微分演算行列



"""DHparamを更新"""
function update_DHparams(q, DHparams=DHparams_neutral)
    #q[2] = q[2] + pi/2
    for i in 1:7
        if i == 2
            DHparams[i].θ = q[i] + pi/2
        else
            DHparams[i].θ = q[i]
        end
    end
    DHparams
end



"""
    HTM(p::DHparam)
同時変換行列(i-1)T(i)
"""
function HTM(p::DHparam{T}) where T
    [
        cos(p.θ) -sin(p.θ) 0.0 p.a
        sin(p.θ)*cos(p.α) cos(p.θ)*cos(p.α) -sin(p.α) -p.d*sin(p.α)
        sin(p.θ)*sin(p.α) cos(p.θ)*sin(p.α) cos(p.α) p.d*cos(p.α)
        0.0 0.0 0.0 1.0
    ]
end


"""同時変換行列を全部計算"""
function calc_HTMs_local_and_global(DHparams::Vector{DHparam{T}}) where T
    HTMs_local = Vector{Matrix{T}}(undef, 10)
    HTMs_global = Vector{Matrix{T}}(undef, 10)
    for i in 1:10
        if i == 1
            HTMs_local[i] = HTM_BL_Wo
            HTMs_global[i] = HTM_BL_Wo
        elseif i == 2
            HTMs_local[i] = HTM_0_BL
            HTMs_global[i] =  HTMs_global[i-1] * HTM_0_BL
        elseif i == 10
            HTMs_local[i] = HTM_GR_7
            HTMs_global[i] = HTMs_global[i-1] * HTM_GR_7
        else
            HTMs_local[i] = HTM(DHparams[i-2])
            HTMs_global[i] = HTMs_global[i-1] * HTMs_local[i]
        end
    end
    HTMs_local, HTMs_global
end


"""微分？同時変換行列を全部計算"""
function calc_dHTMs(
    HTMs_local::Vector{Matrix{T}},
    HTMs_global::Vector{Matrix{T}},
) where T
    _dTdq_all = Vector{Vector{Matrix{T}}}(undef, 10)
    for i in 1:7
        _dTdq  = Vector{Matrix{T}}(undef, 10)
        for j in 1:8
            if j < i
                _dTdq[j] = zeros(T, 4, 4)
            elseif j == i
                _dTdq[j] = HTMs_global[j+2] * HTM_A
            else
                _dTdq[j] = _dTdq[j-1] * HTMs_local[j+2]
            end
        end
        _dTdq_all[i] = _dTdq
    end
    
    Jax_all = Vector{Matrix{T}}(undef, 8)
    Jay_all = Vector{Matrix{T}}(undef, 8)
    Jaz_all = Vector{Matrix{T}}(undef, 8)
    Jo_all = Vector{Matrix{T}}(undef, 8)

    for i in 1:8
        _Jax = Matrix{T}(undef, 4, 7)
        _Jay = Matrix{T}(undef, 4, 7)
        _Jaz = Matrix{T}(undef, 4, 7)
        _Jo = Matrix{T}(undef, 4, 7)
        for j in 1:7
            _Jax[:, j] = _dTdq_all[j][i][:, 1]
            _Jay[:, j] = _dTdq_all[j][i][:, 2]
            _Jaz[:, j] = _dTdq_all[j][i][:, 3]
            _Jo[:, j] = _dTdq_all[j][i][:, 4]
        end
        Jax_all[i] = _Jax
        Jay_all[i] = _Jay
        Jaz_all[i] = _Jaz
        Jo_all[i] = _Jo
    end
    Jax_all, Jay_all, Jaz_all, Jo_all
end


"""各制御点のヤコビ行列を計算"""
function calc_jacobians(
    Jax_all::Vector{Matrix{T}},
    Jay_all::Vector{Matrix{T}},
    Jaz_all::Vector{Matrix{T}},
    Jo_all::Vector{Matrix{T}},
) where T

    """Joのヤコビ行列"""
    function _calc_Jo_global(Jax, Jay, Jaz, Jo, r_bar)
        z_bar = Jax * r_bar[1,1] + Jay * r_bar[2,1] + Jaz * r_bar[3,1] + Jo
        return z_bar[1:3, :]
    end
    
    # ジョイント基底のヤコビ行列を計算
    Jos_joint_all = Vector{Matrix{T}}(undef, 8)
    for i in 1:8
        Jos_joint_all[i] = _calc_Jo_global(
            Jax_all[i], Jay_all[i], Jaz_all[i], Jo_all[i], zeros(T, 4),
        )
    end

    # 制御点の位置に関するヤコビ行列を計算
    Jos_cpoints_all = Vector{Vector{Matrix{T}}}(undef, 8)
    for i in 1:8
        n = length(cpoints_local[i])
        _J = Vector{Matrix{T}}(undef, n)
        for j in 1:n
            _J[j] = _calc_Jo_global(
                Jax_all[i], Jay_all[i], Jaz_all[i], Jo_all[i], cpoints_local[i][j],
            )
        end
        Jos_cpoints_all[i] = _J
    end

    return Jos_joint_all, Jos_cpoints_all
end


"""制御点の位置と速度を計算

"""
function calc_cpoint_x_and_dx_global(
    HTMs_global::Vector{Matrix{T}},
    Jos_cpoints_all::Vector{Vector{Matrix{T}}},
    dq::Vector{T},
) where T
    cpoints_x_global = Vector{Vector{Vector{T}}}(undef, 8)
    for i in 1:8
        n = length(cpoints_local[i])
        _c = Vector{Vector{T}}(undef, n)
        for j in 1:n
            _c[j] = (HTMs_global[i+2] * cpoints_local[i][j])[1:3]
            #println(typeof(_c[j]))
        end
        cpoints_x_global[i] = _c
    end

    cpoints_dx_global = Vector{Vector{Vector{T}}}(undef, 8)
    for i in 1:8
        n = length(cpoints_local[i])
        _dx = Vector{Vector{T}}(undef, n)
        for j in 1:n
            _dx[j] = Jos_cpoints_all[i][j] * dq
        end
        cpoints_dx_global[i] = _dx
    end

    return cpoints_x_global, cpoints_dx_global
end


"""原点+ジョイント位置を取得（図示用）"""
function calc_joints_x_and_dx_global(
    HTMs_global::Vector{Matrix{T}},
    Jos_joint_all::Vector{Matrix{T}},
    dq::Vector{T},
) where T
    joints_x_global = Vector{Vector{T}}(undef, 11)
    joints_x_global[1] = zeros(T, 3)
    for i in 2:length(joints_x_global)
        joints_x_global[i] = HTMs_global[i-1][1:3, 4]
    end

    joints_dx_global = Vector{Vector{T}}(undef, 11)
    joints_dx_global[1] = zeros(T, 3)
    joints_dx_global[2] = zeros(T, 3)
    joints_dx_global[3] = zeros(T, 3)
    for i in 4:length(joints_dx_global)
        joints_dx_global[i] = Jos_joint_all[i-3] * dq
    end

    return joints_x_global, joints_dx_global
end




"""全部計算"""
function calc_all(q=q_neutral, dq=zeros(Float64, 7),)
    DHparams = update_DHparams(q)
    HTMs_local, HTMs_global = calc_HTMs_local_and_global(DHparams)
    Jax_all, Jay_all, Jaz_all, Jo_all = calc_dHTMs(HTMs_local, HTMs_global)
    Jos_joint_all, Jos_cpoint_all = calc_jacobians(Jax_all, Jay_all, Jaz_all, Jo_all)
    cpoints_x_global, cpoints_dx_global = calc_cpoint_x_and_dx_global(
        HTMs_global, Jos_cpoint_all, dq
    )
    joints_x_global, joints_dx_global = calc_joints_x_and_dx_global(
        HTMs_global, Jos_joint_all, dq
    )
    # for i in 1:8
    #     println("i = ", i)
    #     println(Jos_joint_all[i])
    # end

    return (
        HTMs_local, HTMs_global,
        Jax_all, Jay_all, Jaz_all, Jo_all,
        Jos_joint_all, Jos_cpoint_all,
        cpoints_x_global, cpoints_dx_global,
        joints_x_global, joints_dx_global,
    )
end




#@time for i in 1:1; calc_all() end



function draw_arm(q=q_neutral, dq=zeros(Float64, 7), goal=nothing, obs=nothing)

    _, _, _, _, _, _, _, _, cpoints_x_global, _, joints_x_global, _, = calc_all(q, dq)
    fig = plot()


    x, y, z = split_vec_of_arrays(joints_x_global)
    plot!(
        fig,
        x, y, z,
        aspect_ratio = 1,
        marker=:circle,
        markerα = 0.5,
        label = "joints",
    )

    cname = ("1", "2", "3", "4", "5", "6", "7", "GL")
    for (i, cs) in enumerate(cpoints_x_global)
        x, y, z = split_vec_of_arrays(cs)
        scatter!(
            fig,
            x, y, z,
            label = cname[i]
        )
    end

    if !isnothing(goal)
        scatter!([goal.x[1]], [goal.x[2]], [goal.x[3]], markershape=:star6)
    end

    if !isnothing(obs)
        _o = get_x_from_State(obs)
        x, y, z = split_vec_of_arrays(_o)
        scatter!(x, y, z, markershape=:diamond)
    end

    return fig
end


# function draw_arm(nodes, goal=nothing, obs=nothing)
#     fig = plot()

#     nodes
#     x, y, z = split_vec_of_arrays(joints_x_global)
#     plot!(
#         fig,
#         x, y, z,
#         aspect_ratio = 1,
#         marker=:circle,
#         markerα = 0.5,
#         label = "joints",
#     )

#     cname = ("1", "2", "3", "4", "5", "6", "7", "GL")
#     for (i, cs) in enumerate(cpoints_x_global)
#         x, y, z = split_vec_of_arrays(cs)
#         scatter!(
#             fig,
#             x, y, z,
#             label = cname[i]
#         )
#     end

#     if !isnothing(goal)
#         scatter!([goal[1]], [goal[2]], [goal[3]], markershape=:star6)
#     end

#     if !isnothing(obs)
#         x, y, z = split_vec_of_arrays(obs)
#         scatter!(x, y, z, markershape=:diamond)
#     end

#     return fig
# end


@time fig = draw_arm()


# function draw_arm(fig, q, DHparams, name)
#     """アームをplot by Plots
#     name : trueならright、falseならleft
#     """
#     DHparams = update_DHparams(DHparams, q)
#     os = []
#     Ts = []
#     push!(os, [0.0, 0.0, 0.0])

#     if name
#         T_temp = HTM_BR_Wo
#         push!(Ts, T_temp)
#         push!(os, T_temp[1:3, 4])

#         T_temp = T_temp * HTM_0_BR
#         push!(Ts, T_temp)
#         push!(os, T_temp[1:3, 4])
#     else
#         T_temp = HTM_BL_Wo
#         push!(Ts, T_temp)
#         push!(os, T_temp[1:3, 4])

#         T_temp = T_temp * HTM_0_BL
#         push!(Ts, T_temp)
#         push!(os, T_temp[1:3, 4])
#     end

#     for i in 1:7
#         T_temp = T_temp * T(DHparams[i])
#         push!(Ts, T_temp)
#         push!(os, T_temp[1:3, 4])
#     end

#     T_temp = T_temp * HTM_GR_7
#     push!(Ts, T_temp)
#     push!(os, T_temp[1:3, 4])
#     println(T_temp[1:3, 4])

#     x, y, z = split_vec_of_arrays(os)
#     if name
#         s = "R-"
#     else
#         s = "L-"
#     end

#     plot!(
#         fig,
#         x, y, z,
#         aspect_ratio = 1,
#         marker=:circle,
#         markerα = 0.5,
#         label = s * "joints",
#         )

#     cs_global_all = []
#     for (i, cs_local) in enumerate(c_points_all)
#         cs_global = []
#         for r in cs_local
#             o = Ts[i+2] * r
#             push!(cs_global, o[1:3, :])
#         end
#         push!(cs_global_all, cs_global)
#     end
#     cname = (
#         "1", "2", "3", "4", "5", "6", "7", "GL"
#     )

#     for (i, cs) in enumerate(cs_global_all)
#         if i == 8
#             continue
#         else
#             x, y, z = split_vec_of_arrays(cs)
#             scatter!(
#                 fig, 
#                 x, y, z, label = s * cname[i],
#                 aspect_ratio = 1,
#             )
#         end
#     end
#     fig
# end

# ### neutralを図示 ###
# fig = plot(aspect_ratio = 1,)
# @time fig = draw_arm(fig, qr, DHparams_r, true)
# @time fig = draw_arm(fig, ql, DHparams_l, false)

# # # # 回転
# # # @gif for i in range(0, stop = 360*2, length = 100)
# # #     plot!(fig, camera = (i, 0),)
# # # end
