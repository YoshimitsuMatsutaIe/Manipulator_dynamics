"""RMP制御器"""


using LinearAlgebra
eye(T::Type, n) = Diagonal{T}(I, n)
eye(n) = eye(Float64, n)


"""pullback演算"""
function pullback(f, M, J, dJ=nothing, dx=nothing)
    #pulled_f = zero(f)
    #pulled_M = zero(M)
    if isnothing(dJ) && isnothing(dx)
        pulled_f = J' * f
    else
        pulled_f = J' * (f - M * dJ * dx)
    end
    pulled_M = J' * M * J
    
    return pulled_f, pulled_M
end



"""ソフト正規化関数"""
function soft_normal(v, alpha)
    v_norm = norm(v)
    softmax = v_norm + 1/alpha * log(1 + exp(-2 * alpha * v_norm))
    return v / softmax
end

"""空間を一方向に伸ばす計量"""
function metric_stretch(v, alpha)
    xi = soft_normal(v, alpha)
    return xi * xi'
end

"""基本の計量"""
function basic_metric_H(f, alpha, beta)
    return beta * metric_stretch(f, alpha) + (1 - beta) * eye(3)
end

"""アフィン変換されたシグモイド写像"""
sigma_L(q, q_min, q_max) = (q_max - q_min) * (1 / (1 + exp.(-q))) + q_min

"""ジョイント制限に関する対角ヤコビ行列"""
function D_sigma(q, q_min, q_max)
    diags = (q_max - q_min) * (exp.(-q) / (1 + exp.(-q))^.2)
    return diagm(diags)
end



struct OriginalRMPAttractor{T}
    max_speed::T
    gain::T
    ddq_damp_r::T
    simga_W::T
    sigma_H::T
    metric_damp_r::T
end


"""目標加速度 from OirginalRMP"""
function ddq(p::OriginalRMPAttractor{T}, z, dz, z0) where T
    damp = p.gain / p.max_speed
    a = p.gain * soft_normal(z0-z, p.r) - damp*dz
    return a
end

"""目標計量  from OirginalRMP"""
function inertia_matrix(p::OriginalRMPAttractor{T}, z, dz, z0, ddq) where T
    dis = norm(z0 - z)
    weight = exp(-dis / p.simga_W)
    beta = 1 - exp(-1/2 * (dis / p.sigma_H)^2)
    return weight * basic_metric_H(ddq, p.ddq_damp_r, beta)
end

"""canonical form []"""
function get_canonical(p::OriginalRMPAttractor{T}, z, dz, z0) where T
    a = ddq(p, z, dz, z0)
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
function ddq(p::OriginalRMPCollisionAvoidance{T}, z, dz, z0) where T
    
    x = z - z0
    d = norm(x)
    ∇d = x / d

    # 斥力項
    α_rep = p.gain * exp(-d / p.scale_rep)
    ddq_rep = α_rep * ∇dis

    # ダンピング項
    P_obs = max(0, -dz' * ∇d) * ∇d * ∇d' * dz
    damp_gain = p.gain * p.ratio
    α_damp = damp_gain / (d / p.scale_damp + 1e-7)
    ddq_damp = α_damp * P_obs

    return ddq_rep + ddq_damp
end
