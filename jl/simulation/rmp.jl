"""RMP制御器"""


#using Parameters

using LinearAlgebra
using ForwardDiff  # 自動微分パッケージ?


"""pullback演算"""
function pullbacked_rmp(f, M, J, dJ=nothing, dx=nothing)
    #pulled_f = zero(f)
    #pulled_M = zero(M)
    if isnothing(dJ) && isnothing(dx)
        pulled_f = J' * f
    else
        pulled_f = J' * (f .- M * dJ * dx)
    end
    pulled_M = J' * M * J
    
    return pulled_f, pulled_M
end


# mutable struct NaturalFormRMP{T}
#     f::Vector{T}
#     M::Matrix{T}
# end

# mutable struct CanonicalFormRMP{T}
#     a::Vector{T}
#     M::Matrix{T}
# end



#### オリジナルのRMP ###
"""ソフトマックス関数"""
soft_max(s, α) = s + 1/α * log(1 + exp(-2 * α * s))

"""ソフト正規化関数"""
function soft_normal(v, alpha)
    v_norm = norm(v)
    return v ./ soft_max(v_norm, alpha)
end

"""空間を一方向に伸ばす計量"""
function metric_stretch(v, alpha)
    xi = soft_normal(v, alpha)
    return xi * xi'
end

"""基本の計量"""
function basic_metric_H(f::Vector{T}, alpha::T, beta::T) where T
    return beta .* metric_stretch(f, alpha) + (1 - beta) .* Matrix{T}(I, 3, 3)
end

# """アフィン変換されたシグモイド写像"""
# sigma_L(q, q_min, q_max) = (q_max - q_min) * (1 / (1 + exp.(-q))) + q_min

"""ジョイント制限に関する対角ヤコビ行列"""
function D_sigma(q::Vector{T}, q_min::Vector{T}, q_max::Vector{T}) where T
    diags = (q_max .- q_min) .* (exp.(-q) ./ (1 .+ exp.(-q)).^2)
    #println(diags)
    return diagm(diags)
end


"""OriginalRMPの目標到達制御器のパラメータ"""
struct OriginalRMPAttractor{T}
    max_speed::T
    gain::T
    ddq_damp_r::T
    simga_W::T
    sigma_H::T
    metric_damp_r::T
end



"""目標加速度 from OirginalRMP"""
function ddz(p::OriginalRMPAttractor{T}, z::Vector{T}, dz::Vector{T}, z0::Vector{T}) where T
    damp = p.gain / p.max_speed
    a = p.gain .* soft_normal(z0.-z, p.metric_damp_r) .- damp*dz
    return a
end

"""目標計量  from OirginalRMP"""
function inertia_matrix(
    p::OriginalRMPAttractor{T}, z::Vector{T}, dz::Vector{T}, z0::Vector{T}, ddq::Vector{T}
) where T
    dis = norm(z0 .- z)
    weight = exp(-dis ./ p.simga_W)
    beta = 1 - exp(-1/2 * (dis / p.sigma_H)^2)
    return weight .* basic_metric_H(ddq, p.ddq_damp_r, beta)
end

"""canonical form []"""
function get_canonical(p::OriginalRMPAttractor{T}, z, dz, z0) where T
    a = ddz(p, z, dz, z0)
    M = inertia_matrix(p, z, dz, z0, a)
    return a, M
end


"""natural form ()"""
function get_natural(p::OriginalRMPAttractor{T}, z, dz, z0) where T
    a, M = get_canonical(p, z, dz, z0)
    f = M * a
    return f, M
end


struct OriginalRMPCollisionAvoidance{T}
    scale_rep::T
    scale_damp::T
    ratio::T
    gain::T
    r::T
end

"""障害物回避加速度  from OirginalRMP"""
function ddz(p::OriginalRMPCollisionAvoidance{T}, z, dz, z0) where T
    
    x = z .- z0
    d = norm(x)
    ∇d = x ./ d  # 勾配

    # 斥力項
    α_rep = p.gain .* exp(-d / p.scale_rep)
    ddq_rep = α_rep .* ∇d

    # ダンピング項
    P_obs = max(0.0, -dz' * ∇d) * ∇d * ∇d' * dz
    damp_gain = p.gain .* p.ratio
    α_damp = damp_gain ./ (d / p.scale_damp + 1e-7)
    ddq_damp = α_damp * P_obs

    return ddq_rep .+ ddq_damp
end

"""障害物計量 from OriginalRMP"""
function inertia_matrix(p::OriginalRMPCollisionAvoidance{T}, z, dz, z0, ddq) where T
    d = norm(z .- z0)
    weight = (d / p.r)^2 - 2 * d / p.r + 1
    return weight .* Matrix{T}(I, 3, 3)
end

"""canonical form []"""
function get_canonical(p::OriginalRMPCollisionAvoidance{T}, z, dz, z0) where T
    a = ddz(p, z, dz, z0)
    M = inertia_matrix(p, z, dz, z0, a)
    return a, M
end

"""natural form ()"""
function get_natural(p::OriginalRMPCollisionAvoidance{T}, z, dz, z0) where T
    a, M = get_canonical(p, z, dz, z0)
    f = M * a
    return f, M
end


struct OriginalJointLimitAvoidance{T}
    γ_p::T
    γ_d::T
    λ::T
end

"""ジョイント制限回避加速度 from OriginalRMP"""
function ddz(p::OriginalJointLimitAvoidance{T}, q, dq, q_max, q_min) where T
    z = p.γ_p * (-q) - p.γ_d * dq
    a = inv(D_sigma(q, q_min, q_max)) * z
    return a
end

"""ジョイント制限回避計量 from OriginalRMP"""
function inertia_matrix(p::OriginalJointLimitAvoidance{T}, q, dq) where T
    return p.λ .* Matrix{T}(I, 7, 7)
end

"""canonical form []"""
function get_canonical(p::OriginalJointLimitAvoidance{T}, q, dq, q_max, q_min) where T
    a = ddz(p, q, dq, q_max, q_min)
    M = inertia_matrix(p, q, dq)
    return a, M
end


"""natural form ()"""
function get_natural(p::OriginalJointLimitAvoidance{T}, q, dq, q_max, q_min) where T
    a, M = get_canonical(p, q, dq, q_max, q_min)
    # println(a)
    # println(M)
    f = M * a
    return f, M
end




### fromGDS ###


"""fromGDSのアトラクター　パラメータ"""
struct RMPfromGDSAttractor{T}
    max_speed::T
    gain::T
    f_α::T
    σ_α::T
    σ_γ::T
    wᵤ::T
    wₗ::T
    α::T
    ϵ::T
end

"""目標吸引ポテンシャル2"""
potential_2(x, η) = soft_normal(x, η)

"""ポテンシャルの勾配"""
function ∇potential_2(x, η)
    x_norm = norm(x)
    return (1-exp(-2*η*x_norm)) / (1+exp(-2*η*x_norm)) * x / x_norm
end

"""?"""
α_or_γ(x, σ) = exp(-(norm(x))^2 / (2*σ^2))

"""重み行列（fromGDSのアトラクターで使用）"""
w(x, σ_γ, wᵤ, wₗ) = (wᵤ-wₗ) * α_or_γ(x, σ_γ) + wₗ

"""fromGDSのアトラクター慣性行列"""
function inertia_matrix(x, p::RMPfromGDSAttractor{T}, x₀) where T
    z = x .- x₀
    ∇pot = ∇potential_2(z, p.α)
    α = α_or_γ(z, p.σ_α)
    return w(z, p.σ_γ, p.wᵤ, p.wₗ) .* ((1-α) .* ∇pot * ∇pot' .+ (α + p.ϵ).*Matrix{T}(I, 3, 3))
end

"""力用"""
function xMx(x, p::RMPfromGDSAttractor{T}, ẋ, x₀) where T
    ẋ' * inertia_matrix(x, p, x₀) * ẋ
end

"""曲率項"""
function ξ(p::RMPfromGDSAttractor{T}, x, ẋ, x₀) where T
    # 第一項を計算
    A = Matrix{T}(undef, 3, 3)
    for i in 1:3
        _M(x) = inertia_matrix(x, p, x₀)[:, i]
        _jacobian_M = ForwardDiff.jacobian(_M, x)
        #println(_jacobian_M)
        A[:, i] = _jacobian_M * ẋ
    end
    A *= ẋ

    # 第二項を計算
    _xMx(x) = xMx(x, p, ẋ, x₀)
    B = 1/2 .* ForwardDiff.gradient(_xMx, x)  # 便利
    
    return A .- B
end

"""fromGDSのアトラクター力"""
function f(p::RMPfromGDSAttractor{T}, x, ẋ, x₀, M) where T
    z = x .- x₀
    damp = p.gain / p.max_speed
    return M * (-p.gain .* soft_normal(z, p.f_α) .- damp .* ẋ) .- ξ(p, x, ẋ, x₀)
end

""""""
function get_natural(p::RMPfromGDSAttractor{T}, x, ẋ, x₀) where T
    M = inertia_matrix(x, p, x₀)
    return f(p, x, ẋ, x₀, M), M
end


# # テスト
# p = RMPfromGDSAttractor(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0)
# x = [1.0, 2.0, 3.0]
# dx = [1.0, 2.0, 3.0]
# x0 = [0.0, 0.0, 0.0]

# M = inertia_matrix(x, p, x0)
# temp = ξ(p, x, dx, x0, M)


struct RMPfromGDSCollisionAvoidance{T}
    rw::T
    σ::T
    α::T
end

"""重み関数"""
w(s) = s^(-4)

#function w(s)

"""重み関数の微分"""
dwds(s) = -4 * s^(-5)

function u(ṡ, σ)
    if ṡ < 0.0
        return 1 - exp(-ṡ^2 / (2*σ^2))
    else
        return 0.0
    end
end

function dudṡ(ṡ, σ)
    if ṡ < 0.0
        return -exp(-ṡ^2 / (2*σ^2)) * (-ṡ / σ^2)
    else
        return 0.0
    end
end

δ(s, ṡ, σ) = u(ṡ, σ) + 1/2 * ṡ * dudṡ(ṡ, σ)

"""曲率項"""
ξ(s, ṡ, σ) = 1/2 * u(ṡ, σ) * dwds(s) * ṡ^2

"""障害物回避ポテンシャル"""
Φ₁(s, α) = -1/2 * α * s^(-2)

"""障害物回避ポテンシャルの勾配"""
∇Φ₁(s, α) = α * s^(-3)


"""fromGDSの障害物回避力"""
function f(p::RMPfromGDSCollisionAvoidance{T}, s, ṡ) where T
    return -w(s) * ∇Φ₁(s, p.α) - ξ(s, ṡ, p.σ)
end

"""fromGDSの障害物回避慣性行列"""
function inertia_matrix(p::RMPfromGDSCollisionAvoidance{T}, s, ṡ) where T
    return w(s) * δ(s, ṡ, p.σ)
end

function get_natural(p::RMPfromGDSCollisionAvoidance{T}, x, ẋ, x₀, ẋ₀=zeros(Float64, 3)) where T
    s_vec = x₀ .- x
    ṡ_vec = ẋ₀ .- ẋ
    s = norm(s_vec)
    ṡ = (1/s .* dot(s_vec, ṡ_vec))[1]
    println("x = ", x)
    println("dx = ", ẋ)
    println("s = ", s)
    println("s_dot_vec = ", ṡ_vec)

    m = inertia_matrix(p, s, ṡ)
    _f = f(p, s, ṡ)

    J = -(x .- ẋ)' ./ s
    J̇ = -s^(-2) .* (ṡ_vec' .- s_vec' .* ṡ)

    _f, M = pullbacked_rmp(_f, m, J, J̇, ẋ)
    println(_f)
    println(M)
    println()
    return _f, M
end

