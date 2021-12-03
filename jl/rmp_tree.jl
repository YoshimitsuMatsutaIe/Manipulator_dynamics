

# """
# rmp-treeで所望の加速度を計算
# """
# module RMPTree


# export update_nodes
# export calc_desired_ddq

#push!(LOAD_PATH, ".")

using LinearAlgebra



include("utils.jl")

include("rmp.jl")
include("kinematics.jl")


using .Kinematics: calc_all, q_max, q_min, q_neutral, cpoints_local
# using .RMP
#using .Utilis






# """所望の加速度を計算（resolve演算結果を返す）

# # args
# - nodes :  
# - rmp_param :  
# - goal :  
# - obs :  
# """
# function calc_desired_ddq(
#     nodes::Vector{Vector{Node{T}}},
#     rmp_param::NamedTuple,
#     goal::State{T},
#     obs::Vector{State{T}}
#     ) where T

#     root_f = zeros(T, 7)
#     root_M = zeros(T, 7, 7)

#     for i in 1:9
#         if i == 1  # ジョイント制限rmp
#             _f, _M = get_natural(
#                 rmp_param.joint_limit_avoidance, nodes[i][1].x, nodes[i][1].dx,
#                 q_max, q_min, q_neutral
#             )
#             #_f = zeros(T, 7)
#             #_M = zeros(T, 7, 7)
#             @. root_f += _f
#             @. root_M += _M
            
#         else
#             if i == 9  # 目標王達rmpを追加
#             _f, _M = get_natural(
#                 rmp_param.attractor, nodes[i][1].x, nodes[i][1].dx, goal.x
#             )
#             _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][1].Jo,)
#             @. root_f += _pulled_f
#             @. root_M += _pulled_M
#             end

#             # 障害物回避rmpを追加
#             for j in 1:length(nodes[i])
#                 for k in 1:length(obs)
#                     _f, _M = get_natural(
#                         rmp_param.obs_avoidance[i-1], nodes[i][j].x, nodes[i][j].dx, obs[k].x
#                     )
#                     _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][j].Jo,)
#                     @. root_f += _pulled_f
#                     @. root_M += _pulled_M
#                 end
#             end
#         end
#     end

#     #root_f += zeros(T, 7)
#     #root_M += zeros(T, 7, 7)

#     #println("root_f = ", root_f)
#     #println("root_M = ", root_M)


#     ddq = pinv(root_M) * root_f
#     #println("ddq = ", ddq)
#     #ddq = np.linalg.pinv(root_M) * root_f

#     #ddq = zeros(T, 7)
#     return ddq
# end


"""nodeを新規作成"""
function update_nodes(nodes::Nothing, q::Vector{T}, dq::Vector{T}) where T
    # 更新
    _, _,
    _, _, _, _,
    _, Jos_cpoint_all,
    cpoints_x_global, cpoints_dx_global,
    joints_x_global, joints_dx_global = calc_all(q, dq)
    nodes = Vector{Vector{Node{T}}}(undef, 9)
    for i in 1:9
        #println("i = ", i)
        if i == 1  # rootノード
            _children_node = Vector{Node{T}}(undef, 1)
            _children_node[1] = Node(
                x = q,
                dx = dq,
                Jo = Matrix{T}(I, 7, 7),
                f = zeros(T, 7),
                M = zeros(T, 7, 7)
                )
            #println("wow")
        else
            n = length(cpoints_local[i-1]) + 1
            _children_node = Vector{Node{T}}(undef, n)
            _children_node[1] = Node(
                x = joints_x_global[i+2],
                dx = joints_dx_global[i+2],
                Jo = nothing,
                f = nothing,
                M = nothing
            )  # ジョイント点

            for j in 2:n  # 制御点
                _children_node[j] = Node(
                    x = cpoints_x_global[i-1][j-1],
                    dx = cpoints_dx_global[i-1][j-1],
                    Jo = Jos_cpoint_all[i-1][j-1],
                    f = zeros(T, 3),
                    M = zeros(T, 3, 3),
                )
            end
        end
        nodes[i] = _children_node
    end
    return nodes
end



"""nodeを更新"""
function update_nodes!(;
    nodes::Vector{Vector{Node{T}}}, q::Vector{T}, dq::Vector{T},
    rmp_param::NamedTuple, goal::State{T}, obs::Vector{State{T}}
    ) where T
    # 更新
    _, _,
    _, _, _, _,
    _, Jos_cpoint_all,
    cpoints_x_global, cpoints_dx_global,
    joints_x_global, joints_dx_global = calc_all(q, dq)

    #hoge = nodes[1][1].x
    #nodes = copy(nodes)
    #println(old[1][1].x)

    for i in 1:9
        for j in 1:length(nodes[i])
            if i == 1  # root
                nodes[i][j].x = q
                nodes[i][j].dx = dq
                nodes[i][j].Jo = Matrix{T}(I, 7, 7)
                nodes[i][j].f, nodes[i][j].M = get_natural(
                    rmp_param.joint_limit_avoidance, nodes[i][1].x, nodes[i][1].dx,
                    q_max, q_min, q_neutral
                )

                # # ジョイント制限を考慮しない
                # nodes[i][j].f = zeros(T, 7)
                # nodes[i][j].M = zeros(T, 7, 7)


            else  # root以外
                if j == 1   # ジョイント位置
                    nodes[i][j].x = joints_x_global[i+2]
                    nodes[i][j].dx = joints_dx_global[i+2]
                    nodes[i][j].Jo = nothing
                    nodes[i][j].f = nothing
                    nodes[i][j].M = nothing
                else
                    nodes[i][j].x = cpoints_x_global[i-1][j-1]
                    nodes[i][j].dx = cpoints_dx_global[i-1][j-1]
                    nodes[i][j].Jo = Jos_cpoint_all[i-1][j-1]

                    f = zeros(T, 3)
                    M = zeros(T, 3, 3)
                    for k in 1:length(obs)
                        _f, _M = get_natural(
                            rmp_param.obs_avoidance[i-1], nodes[i][j].x, nodes[i][j].dx, obs[k].x
                        )
                        @. f += _f
                        @. M += _M
                    end
                    nodes[i][j].f = f
                    nodes[i][j].M = M


                    if i==9  # エンドエフェクタのときはさらに目標rmpを追加
                        _f, _M = get_natural(
                            rmp_param.attractor, nodes[i][j].x, nodes[i][j].dx, goal.x
                        )
                        @. nodes[i][j].f += _f
                        @. nodes[i][j].M += _M
                    end
                end
            end
        end
    end
    # println(nodes[1][1].x)
    # println(nodes[1][1].x == old[1][1].x)
    return nodes
end



"""nodesからf, Mを集めて所望の加速度を生成

新手法  
"""
function calc_desired_ddq(nodes::Vector{Vector{Node{T}}}) where T
    root_f = zeros(T, 7)
    root_M = zeros(T, 7, 7)

    for i in 1:9
        for j in 1:length(nodes[i])  # 目標rmpと障害物rmpの合算
            if i == 1 && !isnothing(nodes[i][j].f)
                _pulled_f, _pulled_M = pullbacked_rmp(
                    nodes[i][j].f, nodes[i][j].M, nodes[i][j].Jo,
                )
                @. root_f += _pulled_f
                @. root_M += _pulled_M
            elseif i > 1 && !isnothing(nodes[i][j].f)
                _pulled_f, _pulled_M = pullbacked_rmp(
                    nodes[i][j].f, nodes[i][j].M, nodes[i][j].Jo,
                )
                @. root_f += _pulled_f
                @. root_M += _pulled_M
            end
        end
    end

    #root_f += zeros(T, 7)
    #root_M += zeros(T, 7, 7)

    #println("root_f = ", root_f)
    #println("root_M = ", root_M)
    ddq = pinv(root_M) * root_f
    #println("ddq = ", ddq)
    #ddq = np.linalg.pinv(root_M) * root_f

    #ddq = zeros(T, 7)
    return ddq
end




"""nodesからジョイント原点位置を取得（アニメ用）"""
function get_joint_positions_from_nodes(nodes::Vector{Vector{Node{T}}}) where T

    xs = Vector{Vector{T}}(undef, 8)

    for i in 2:9
        xs[i-1] = nodes[i][1].x
    end
        return xs
end


"""nodesから制御点位置を取得（アニメ用）"""
function get_cpoints_from_nodes(nodes::Vector{Vector{Node{T}}}) where T

    xs = Vector{Vector{Vector{T}}}(undef, 8)
    for i in 2:9
        _xs = Vector{Vector{T}}(undef, length(nodes[i])-1)
        for j in 2:length(nodes[i])
            _xs[j-1] = nodes[i][j].x
        end
        xs[i-1] = _xs
    end

    return xs

end


"""nodesからrmpを取得（グラフアニメ用）"""
function get_rmps_from_nodes(nodes::Vector{Vector{Node{T}}}) where T
    fs = Vector{Vector{Vector{T}}}(undef, 9)
    Ms = Vector{Vector{Matrix{T}}}(undef, 9)


    # fs[1][1] = nodes[1][1].f
    # Ms[1][1] = nodes[1][1].M

    for i in 1:9
        if i == 1
            _fs = Vector{Vector{T}}(undef, length(nodes[i]))
            _Ms = Vector{Matrix{T}}(undef, length(nodes[i]))
            _fs[i] = nodes[1][1].f
            _Ms[i] = nodes[1][1].M
            
        else
            
            _fs = Vector{Vector{T}}(undef, length(nodes[i]) - 1)
            _Ms = Vector{Matrix{T}}(undef, length(nodes[i]) - 1)
            for j in 2:length(nodes[i])
                _fs[j-1] = nodes[i][j].f
                _Ms[j-1] = nodes[i][j].M
            end
        end
        fs[i] = _fs
        Ms[i] = _Ms
    end

    return fs, Ms
end

"""fsからfのノルムだけ取得"""
function get_f_norms_from_fs(fs::Vector{Vector{Vector{T}}}) where T
    f_norms = Vector{Vector{T}}(undef, 9)

    for (i, _fs) in enumerate(fs)
        _f_norms = Vector{T}(undef, length(_fs))
        for (j, f) in enumerate(_fs)
            _f_norms[j] = norm(f)
        end
        f_norms[i] = _f_norms
    end
    return f_norms
end


"""fノルムsを時刻でまとめる"""
function get_f_norms_from_nodes_list(
    t_max::Int64,
    nodes_list::Vector{Vector{Vector{Node{T}}}}
    ) where T

    f_norms_list = Vector{Vector{Vector{T}}}(undef, 9)

    for (t, nodes) in enumerate(nodes_list)
        fs, _ = get_rmps_from_nodes(nodes)
        f_norms = get_f_norms_from_fs(fs)

        if t == 1  # 1回目は配列を宣言する
            
            for (i, f) in enumerate(f_norms)
                _f_norms_list = Vector{Vector{T}}(undef, length(f))
                for (j, _f) in enumerate(f)
                    _f_norms_list[j] = Vector{T}(undef, t_max)
                    _f_norms_list[j][1] = _f
                end
                f_norms_list[i] = _f_norms_list
            end
        else
            for (i, f) in enumerate(f_norms)
                for (j, _f) in enumerate(f)
                    f_norms_list[i][j][t] = _f
                end
            end
        end
    end
    return f_norms_list
end

#end