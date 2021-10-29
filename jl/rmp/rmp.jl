"""RMP制御器"""


# """
# RMPいろいろ
# """
# module RMP

using LinearAlgebra
using ForwardDiff  # 自動微分パッケージ


# export pullbacked_rmp
# export  get_natural

# export OriginalRMPAttractor
# export OriginalRMPCollisionAvoidance
# export OriginalJointLimitAvoidance
# export RMPfromGDSAttractor
# export RMPfromGDSCollisionAvoidance


"""pullback演算"""
function pullbacked_rmp(f, M, J, dJ, dx)
    pulled_f = J' * (f .- M * dJ * dx)
    pulled_M = J' * M * J
    return pulled_f, pulled_M
end

"""pullback演算"""
function pullbacked_rmp(f, M, J)
    pulled_f = J' * f
    pulled_M = J' * M * J
    return pulled_f, pulled_M
end




#### オリジナルのRMP ###
"""ソフトマックス関数"""
function soft_max(s::T, α::T) where T
    s + 1/α * log(1 + exp(-2 * α * s))
end

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
    diags = (q_max .- q_min) .* (exp.(-q) ./ (1.0 .+ exp.(-q)).^2)
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
    beta = 1.0 - exp(-1/2 * (dis / p.sigma_H)^2)
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
    P_obs = max(0.0, dot(-dz, ∇d)) * ∇d * ∇d' * dz
    damp_gain = p.gain * p.ratio
    α_damp = damp_gain / (d / p.scale_damp + 1e-7)
    ddq_damp = α_damp .* P_obs

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
    f = M * a * 0.01
    return f, M
end


struct OriginalJointLimitAvoidance{T}
    γ_p::T
    γ_d::T
    λ::T
end

"""ジョイント制限回避加速度 from OriginalRMP"""
function ddz(
    p::OriginalJointLimitAvoidance{T},
    q::Vector{T}, dq::Vector{T}, q_max::Vector{T}, q_min::Vector{T}
) where T
    z = p.γ_p .* (-q) .- p.γ_d .* dq
    a = inv(D_sigma(q, q_min, q_max)) * z
    return a
end

"""ジョイント制限回避計量 from OriginalRMP"""
function inertia_matrix(p::OriginalJointLimitAvoidance{T}, q, dq) where T
    return p.λ .* Matrix{T}(I, 7, 7)
end

"""canonical form []"""
function get_canonical(p::OriginalJointLimitAvoidance{T}, q, dq, q_max, q_min) where T
    return ddz(p, q, dq, q_max, q_min), inertia_matrix(p, q, dq)
end


"""natural form ()"""
function get_natural(p::OriginalJointLimitAvoidance{T}, q, dq, q_max, q_min) where T
    a, M = get_canonical(p, q, dq, q_max, q_min)
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
potential_2 = soft_normal

"""ポテンシャルの勾配"""
function ∇potential_2(x, η)
    x_norm = norm(x)
    return (1-exp(-2*η*x_norm)) / (1+exp(-2*η*x_norm)) .* x ./ x_norm
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
function xMx(x, p::RMPfromGDSAttractor{T}, dx, x₀) where T
    dx' * inertia_matrix(x, p, x₀) * dx
end

"""曲率項"""
function ξ(p::RMPfromGDSAttractor{T}, x, dx, x₀) where T
    # 第一項を計算
    A = Matrix{T}(undef, 3, 3)
    for i in 1:3
        _M(x) = inertia_matrix(x, p, x₀)[:, i]
        _jacobian_M = ForwardDiff.jacobian(_M, x)
        #println(_jacobian_M)
        A[:, i] = _jacobian_M * dx
    end
    A *= dx

    # 第二項を計算
    _xMx(x) = xMx(x, p, dx, x₀)
    B = 1/2 .* ForwardDiff.gradient(_xMx, x)  # 便利
    
    return A .- B
end

"""fromGDSのアトラクター力"""
function f(p::RMPfromGDSAttractor{T}, x, dx, x₀, M) where T
    z = x .- x₀
    damp = p.gain / p.max_speed
    return M * (-p.gain .* soft_normal(z, p.f_α) .- damp .* dx) .- ξ(p, x, dx, x₀)
end

""""""
function get_natural(p::RMPfromGDSAttractor{T}, x, dx, x₀) where T
    M = inertia_matrix(x, p, x₀)
    return f(p, x, dx, x₀, M), M
end


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

function w2(s, rw=1.0)
    if rw - s > 0.0
        return (rw - s)^2 / s
    else
        return 0.0
    end
end

function dw2(s, rw)
    if rw - s > 0.0
        return (-2*(rw - s) * s + (rw - s)) / s^2
    else
        return 0.0
    end
end


function u(ds, σ)
    if ds < 0.0
        return 1 - exp(-ds^2 / (2*σ^2))
    else
        return 0.0
    end
end

function du(ds, σ)
    if ds < 0.0
        return -exp(-ds^2 / (2*σ^2)) * (-ds / σ^2)
    else
        return 0.0
    end
end

δ(s, ds, σ) = u(ds, σ) + 1/2 * ds * du(ds, σ)

"""曲率項"""
ξ(s, ds, σ, rw) = 1/2 * u(ds, σ) * dw2(s, rw) * ds^2

"""障害物回避ポテンシャル"""
Φ1(s, α, rw) = 1/2 * α * w2(s, rw)^2

"""障害物回避ポテンシャルの勾配"""
∇Φ1(s, α, rw) = α * w2(s, rw) * dw2(s, rw)


"""fromGDSの障害物回避力"""
function f(p::RMPfromGDSCollisionAvoidance{T}, s, ds) where T
    return -w2(s, p.rw) * ∇Φ1(s, p.α, p.rw) - ξ(s, ds, p.σ, p.rw)
end

"""fromGDSの障害物回避慣性行列"""
function inertia_matrix(p::RMPfromGDSCollisionAvoidance{T}, s, ds) where T
    return w2(s, p.rw) * δ(s, ds, p.σ)
end

function get_natural(p::RMPfromGDSCollisionAvoidance{T}, x, dx, x₀, dx₀=zeros(Float64, 3)) where T
    s_vec = x₀ .- x
    ds_vec = dx₀ .- dx
    s = norm(s_vec)
    ds = (1/s .* dot(s_vec, ds_vec))[1]


    m = inertia_matrix(p, s, ds)
    _f = f(p, s, ds)

    J = -(x .- dx)' ./ s
    J̇ = -s^(-2) .* (ds_vec' .- s_vec' .* ds)

    _f, M = pullbacked_rmp(_f, m, J, J̇, dx)

    # println("x = ", x)
    # println("dx = ", dx)
    # println("s = ", s)
    # println("s_dot_vec = ", ds_vec)
    # println(_f)
    # println(M)
    # println()
    return _f, M
end

#end