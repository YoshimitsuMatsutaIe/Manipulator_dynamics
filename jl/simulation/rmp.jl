"""RMP制御器"""


#using Parameters

using LinearAlgebra


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

"""ソフト正規化関数"""
function soft_normal(v, alpha)
    v_norm = norm(v)
    softmax = v_norm + 1/alpha * log(1 + exp(-2 * alpha * v_norm))
    return v ./ softmax
end

"""空間を一方向に伸ばす計量"""
function metric_stretch(v, alpha)
    xi = soft_normal(v, alpha)
    return xi * xi'
end

"""基本の計量"""
function basic_metric_H(f, alpha::T, beta::T) where T
    return beta .* metric_stretch(f, alpha) + (1 - beta) .* Matrix{T}(I, 3, 3)
end

# """アフィン変換されたシグモイド写像"""
# sigma_L(q, q_min, q_max) = (q_max - q_min) * (1 / (1 + exp.(-q))) + q_min

"""ジョイント制限に関する対角ヤコビ行列"""
function D_sigma(q, q_min, q_max)
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
function ddz(p::OriginalRMPAttractor{T}, z, dz, z0) where T
    damp = p.gain / p.max_speed
    a = p.gain .* soft_normal(z0.-z, p.metric_damp_r) .- damp*dz
    return a
end

"""目標計量  from OirginalRMP"""
function inertia_matrix(p::OriginalRMPAttractor{T}, z, dz, z0, ddq) where T
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

struct RMPfromGDSAttractor{T}
    max_speed::T
    gain::T
    f_α::T
    σ_α::T
    σ_γ::T
    wᵤ::T
    wₗ::T
    ddq_damp_r::T
    α::T
    ϵ::T
end


struct RMPfromGDSCollisionAvoidance{T}
    rw::T
    σ::T
    α::T
end

"""重み関数"""
w(s) = s^(-4)

"""重み関数の微分"""
dwds(s) = -4 * s^(-5)

function u(ṡ, σ)
    if ṡ < 0.0
        return -exp(-ṡ^2 / (2*σ^2))
    else
        return 0.0
    end
end

function dudṡ(ṡ, σ)
    if ds < 0.0
        return -exp(-ṡ^2 / (2*σ^2)) * (-ṡ / σ^2)
    else
        return 0.0
    end
end

δ(s, ṡ, σ) = u(ṡ, σ) + 1/2 * ṡ * dudṡ(ṡ, σ)

ξ(s, ṡ, σ) = 1/2 * u(ṡ, σ) * dwds(s) * ṡ^2

Φ₁(s, α) = -1/2 * α * s^(-2)

∇Φ₁(s, α) = α * s^(-3)


"""障害物回避力"""
function f(p::RMPfromGDSCollisionAvoidance{T}, s, ṡ) where T
    return -w(s) * ∇Φ₁(s, p.α) - ξ(s, ṡ, p.σ)
end

"""障害物回避慣性行列"""
function inertia_matrix(p::RMPfromGDSCollisionAvoidance{T}, s, ṡ) where T
    return w(s) * δ(s, ṡ, p.σ)
end

function get_natural(p::RMPfromGDSCollisionAvoidance{T}, x, ẋ, x₀, ẋ₀) where T
    s_vec = x₀ .- x
    ṡ_vec = ẋ₀ - ẋ
    s = norm(s_vec)
    ṡ = 1/s * s_vec * (ṡ_vec)

    m = inertia_matrix(p, s, ṡ)
    f = f(p, s, ṡ)

    J = -(x .- ẋ)' ./ s
    J̇ = -s^(-2) .* (ṡ_vec' .- s_vec' .* ṡ)

    f, M = pullbacked_rmp(f, m, J, J̇, ẋ)
    return f, M
end