using CPUTime
using Plots
using LinearAlgebra

include("../utils.jl")

include("./rmp.jl")
include("../kinematics/old.jl")


# function jisaku_solve_euler(dx, x₀, t_span, Δt)
#     """オイラー法
#     ・参考にしました: https://twitter.com/genkuroki/status/1301832571131633665/photo/1
#     """

#     t = range(t_span..., step = Δt)  # 時間軸
#     x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

#     x[1] = x₀  # 初期値
#     for i in 1:length(x)-1
#         x[i+1] = x[i] + dx(t[i], x[i])Δt
#     end

#     t, x
# end


# function jisaku_solve_RungeKutta(dx, x₀, t_span, Δt)
#     """ルンゲクッタ法（4次）"""

#     t = range(t_span..., step = Δt)  # 時間軸
#     x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

#     x[1] = x₀  # 初期値
#     for i in 1:length(x)-1
#         k₁ = dx(t[i], x[i])
#         k₂ = dx(t[i]+Δt/2, x[i]+k₁*Δt/2)
#         k₃ = dx(t[i]+Δt/2, x[i]+k₂*Δt/2)
#         k₄ = dx(t[i]+Δt, x[i]+k₃*Δt)
#         x[i+1] = x[i] + (k₁ + 2k₂ + 2k₃ +k₄)Δt/6
#     end

#     t, x
# end



"""加速度指令を計算（resolve演算結果を返す）"""
function calc_ddq()





    return ddq
end






"""ひとまずシミュレーションやってみｓる"""
function run_simulation()

    const TIME_SPAN = 10
    const Δt = 0.01

    goal = [1; 1; 1]
    obs = [2; 2; 2]





end