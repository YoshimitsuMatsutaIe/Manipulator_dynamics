"""RMP制御器"""


using LinearAlgebra


"""
pullback演算
"""
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


mutable struct OriginalRMPAttractorParam{T}
    max_speed::T
    gain::T
    ddq_damp_r::T
    simga_W::T
    sigma_H::T
    metric_damp_r::T
end

mutable struct OriginalRMPCollisionAvoidanceParam{T}
    scale_rep::T
    scale_damp::T
    ratio::T
    gain::T
    r::T
end

"""ソフト正規化関数"""
function soft_normal(v, alpha)
    v_norm = norm(v)
    softmax = v_norm + 1/alpha * log(1 + exp(-2 * alpha * v_norm))
    return v / softmax
end

function metric_stretch(v, alpha)
    xi = soft_normal(v, alpha)
    return xi * xi'
end

function basic_metric_H(f, alpha, beta)
    return beta * metric_stretch(f, alpha) + (1 - beta) * eye()