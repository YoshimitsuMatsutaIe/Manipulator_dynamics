using CPUTime
using Plots
using LinearAlgebra
using StaticArrays
using ArraysOfArrays

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

const q_neutral = [0; -31; 0; 43; 0; 72; 0] * pi/180  # ニュートラルの姿勢

qr = q_neutral  # 右手の関節角度ベクトル
ql = q_neutral  # 左手の関節角度ベクトル

DHparams_r = [
    DHparam(0.0, 0.0, 0.0, qr[1])
    DHparam(-pi/2, L1, 0.0, qr[2]+pi/2)
    DHparam(pi/2, 0.0, L2, qr[3])
    DHparam(-pi/2, L3, 0.0, qr[4])
    DHparam(pi/2, 0.0, L4, qr[5])
    DHparam(-pi/2, L5, 0.0, qr[6])
    DHparam(pi/2, 0.0, 0.0, qr[7])
]

DHparams_l = [
    DHparam(0.0, 0.0, 0.0, ql[1])
    DHparam(-pi/2, L1, 0.0, ql[2]+pi/2)
    DHparam(pi/2, 0.0, L2, ql[3])
    DHparam(-pi/2, L3, 0.0, ql[4])
    DHparam(pi/2, 0.0, L4, ql[5])
    DHparam(-pi/2, L5, 0.0, ql[6])
    DHparam(pi/2, 0.0, 0.0, ql[7])
]

# 制御点
const c_points_1 = (
    [0.0; L1/2; -L0/2; 1.0],
    [0.0; -L1/2; -L0/2; 1.0],
    [L1/2; 0.0; -L0/2; 1.0],
    [-L1/2; 0.0; -L0/2; 1.0],
)


const c_points_2 = (
    [0.0; 0.0; L3/2; 1.0],
    [0.0; 0.0; -L3/2; 1.0],
)

const c_points_3 = (
    [0.0; L3/2; -L2*2/3; 1.0],
    [0.0; -L3/2; -L2*2/3; 1.0],
    [L3/2; 0.0; -L2*2/3; 1.0],
    [-L3/2; 0.0; -L2*2/3; 1.0],
    [0.0; L3/2; -L2*1/3; 1.0],
    [0.0; -L3/2; -L2*1/3; 1.0],
    [L3/2; 0.0; -L2*1/3; 1.0],
    [-L3/2; 0.0; -L2*1/3; 1.0],
)

const c_points_4 = (
    [0.0; 0.0; L3/2; 1.0],
    [0.0; 0.0; -L3/2; 1.0],
)

const c_points_5 = (
    [0.0; L5/2; -L4/2; 1.0],
    [0.0; -L5/2; -L4/2; 1.0],
    [L5/2; 0.0; -L4/2; 1.0],
    [-L5/2; 0.0; -L4/2; 1.0],
)

const c_points_6 = (
    [0.0; 0.0; L5/2; 1.0],
    [0.0; 0.0; -L5/2; 1.0]
)

const c_points_7 = (
    [0.0; L5/2; L6/2; 1.0],
    [0.0; -L5/2; L6/2; 1.0],
    [L5/2; 0.0; L6/2; 1.0],
    [-L5/2; 0.0; L6/2; 1.0],
)

const c_points_GL = (
    [0.0; 0.0; L6/3; 1.0]
)

const c_points_all = (
    c_points_1, c_points_2, c_points_3, c_points_4, c_points_5, c_points_6, c_points_7, c_points_GL
)


const N_of_points_all = (
    length(c_points_1),
    length(c_points_2),
    length(c_points_3),
    length(c_points_4),
    length(c_points_5),
    length(c_points_6),
    length(c_points_7),
    length(c_points_GL),
)

const N_of_cpoints = sum(N_of_points_all)

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

"""DHparamを更新"""
function update_DHparams(DHparams, q)
    q[2] = q[2] + pi/2
    for i in 1:7
        DHparams[i].θ = q[i]
    end
    DHparams
end

function split_vec_of_arrays(u)
    vec.(u) |>
    x -> VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end


"""同時変換行列(i-1)T(i)"""
function HTM(p::DHparam{T}) where T
    [
        cos(p.θ) -sin(p.θ) 0.0 p.a
        sin(p.θ)*cos(p.α) cos(p.θ)*cos(p.α) -sin(p.α) -p.d*sin(p.α)
        sin(p.θ)*sin(p.α) cos(p.θ)*sin(p.α) cos(p.α) p.d*cos(p.α)
        0.0 0.0 0.0 1.0
    ]
end



function get_HTMs(DHparams::Vector{DHparam{T}}) where T
    HTMs = Vector{Matrix{T}}(undef, 10)
    HTMs_Wo = Vector{Matrix{T}}(undef, 10)
    for i in 1:10
        if i == 1
            HTMs[i] = HTM_BL_Wo
            HTMs_Wo[i] = HTM_BL_Wo
        elseif i == 2
            HTMs[i] = HTM_0_BL
            HTMs_Wo[i] =  HTMs_Wo[i-1] * HTM_0_BL
        elseif i == 10
            HTMs[i] = HTM_GR_7
            HTMs_Wo[i] = HTMs_Wo[i-1] * HTM_GR_7
        else
            HTMs[i] = HTM(DHparams[i-2])
            HTMs_Wo[i] = HTMs_Wo[i-1] * HTMs[i]
        end
    end
    HTMs, HTMs_Wo
end



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
