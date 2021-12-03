"""データ可視化"""



using Plots
#using LaTeXStrings


"""基本的なものだけplot"""
function plot_simulation_data(data, savepath)
    println("基本グラフを作成中...")
    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.q)
    fig_q = plot(data.t, q1, ylabel="q", label="_q1", legend=:outerright)
    plot!(fig_q, data.t, q2, label="_q2")
    plot!(fig_q, data.t, q3, label="_q3")
    plot!(fig_q, data.t, q4, label="_q4")
    plot!(fig_q, data.t, q5, label="_q5")
    plot!(fig_q, data.t, q6, label="_q6")
    plot!(fig_q, data.t, q7, label="_q7")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.dq)
    fig_dq = plot(data.t, q1, ylabel="dq", label="_w1", legend=:outerright)
    plot!(fig_dq, data.t, q2, label="_w2")
    plot!(fig_dq, data.t, q3, label="_w3")
    plot!(fig_dq, data.t, q4, label="_w4")
    plot!(fig_dq, data.t, q5, label="_w5")
    plot!(fig_dq, data.t, q6, label="_w6")
    plot!(fig_dq, data.t, q7, label="_w7")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.ddq)
    fig_ddq = plot(data.t, q1, ylabel="ddq", label="_a1", legend=:outerright)
    plot!(fig_ddq, data.t, q2, label="_a2")
    plot!(fig_ddq, data.t, q3, label="_a3")
    plot!(fig_ddq, data.t, q4, label="_a4")
    plot!(fig_ddq, data.t, q5, label="_a5")
    plot!(fig_ddq, data.t, q6, label="_a6")
    plot!(fig_ddq, data.t, q7, label="_a7")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.desired_ddq)
    fig_desired_ddq = plot(data.t, q1, ylabel="desired ddq", label="a*1", legend=:outerright)
    plot!(fig_desired_ddq, data.t, q2, label="a*2")
    plot!(fig_desired_ddq, data.t, q3, label="a*3")
    plot!(fig_desired_ddq, data.t, q4, label="a*4")
    plot!(fig_desired_ddq, data.t, q5, label="a*5")
    plot!(fig_desired_ddq, data.t, q6, label="a*6")
    plot!(fig_desired_ddq, data.t, q7, label="a*7")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.u)
    fig_u = plot(data.t, q1, ylabel="u", label="u*1", legend=:outerright)
    plot!(fig_u, data.t, q2, label="u*2")
    plot!(fig_u, data.t, q3, label="u*3")
    plot!(fig_u, data.t, q4, label="u*4")
    plot!(fig_u, data.t, q5, label="u*5")
    plot!(fig_u, data.t, q6, label="u*6")
    plot!(fig_u, data.t, q7, label="u*7")


    fig_error = plot(
        data.t, data.error,
        label="err", ylabel="error [m]", legend=:outerright,
        ylim=(0, maximum(data.error))
    )
    fig_dis_to_obs = plot(
        data.t, data.dis_to_obs,
        label="obs", ylabel="min distance to obs [m]", ylims=(0.0,), legend=:outerright
    )


    fig_integrate = plot(
        fig_q, fig_dq, fig_ddq, fig_desired_ddq, fig_u, fig_error, fig_dis_to_obs,
        layout=(7,1),
        size=(500,1400)
    )

    fname = savepath * "basic_his.png"
    savefig(fig_integrate, fname)

    println("基本グラフ作成完了")
end


"""rmpのグラフに使用"""
function styling(i::Int64, j::Int64)
    if i == 1
        n = "ro_"
        ls = :dash
        lw = 5
    elseif i == 9
        n = "ee_"
        ls = :dashdot
        lw = 5
    else
        n = "c" * string(i-1) * "_"
        ls = :solid
        lw = :auto
    end
    n *= string(j)

    return (label=n, linestyle=ls, linewidth=lw)
end

"""rmpの力fをグラフ化"""
function plot_rmp_f(data)
    f_norms_list = get_f_norms_from_nodes_list(
        length(data.q), data.nodes
    )


    fig = plot()
    for (i, fs) in enumerate(f_norms_list)
        for (j, f) in enumerate(fs)
            n = styling(i,j)
            plot!(
                fig,
                data.t, f, label = n.label,
                legend=:outerright, #linestyle=n.linestyle,
                linewidth=n.linewidth
            )
        end
    end

    plot!(
        fig,
        xlabel="time [sec]", ylabel="norm(f)",
        xlims=(0, data.t[end]), ylims=(0,)
    )

    return fig
end


"""rmp関連だけplot"""
function plot_rmp(data, savepath)
    println("rmpグラフを作成中...")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.q)
    fig_q = plot(data.t, q1, ylabel="q", label="_q1_", legend=:outerright)
    plot!(fig_q, data.t, q2, label="_q2_")
    plot!(fig_q, data.t, q3, label="_q3_")
    plot!(fig_q, data.t, q4, label="_q4_")
    plot!(fig_q, data.t, q5, label="_q5_")
    plot!(fig_q, data.t, q6, label="_q6_")
    plot!(fig_q, data.t, q7, label="_q7_")
    plot!(xlims=(0, data.t[end]),)

    fig_error = plot(
        data.t, data.error,
        label="err_", ylabel="error [m]",
        xlims=(0.0, data.t[end]), legend=:outerright,
        ylim=(0, maximum(data.error))
    )
    fig_dis_to_obs = plot(
        data.t, data.dis_to_obs,
        label="obs_", ylabel="min distance to obs [m]",
        xlims=(0, data.t[end]), legend=:outerright
    )

    fig_f = plot_rmp_f(data)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.desired_ddq)
    fig_desired_ddq = plot(data.t, q1, ylabel="desired ddq", label="_a*1", legend=:outerright)
    plot!(fig_desired_ddq, data.t, q2, label="_a*2")
    plot!(fig_desired_ddq, data.t, q3, label="_a*3")
    plot!(fig_desired_ddq, data.t, q4, label="_a*4")
    plot!(fig_desired_ddq, data.t, q5, label="_a*5")
    plot!(fig_desired_ddq, data.t, q6, label="_a*6")
    plot!(fig_desired_ddq, data.t, q7, label="_a*7")

    fig_integrate = plot(
        fig_q, fig_error, fig_dis_to_obs, fig_f, fig_desired_ddq,
        layout=(5,1),
        size=(700,1800)
    )

    fname = savepath * "rmp_his.png"
    savefig(fig_integrate, fname)

    println("rmpグラフを作成完了")
end


"""ロボットアームを描画（q, dqの時刻歴から）"""
function draw_arm(q=q_neutral, dq=zeros(Float64, 7), goal=nothing, obs=nothing, t=nothing, jl=nothing)
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
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]"
    )

    cname = ("1", "2", "3", "4", "5", "6", "7", "GL")
    for (i, cs) in enumerate(cpoints_x_global)
        x, y, z = split_vec_of_arrays(cs)
        scatter!(
            fig,
            x, y, z,
            label = cname[i],
            markersize=5,
        )
    end

    if !isnothing(goal)
        scatter!(
            [goal.x[1]], [goal.x[2]], [goal.x[3]],
            markershape=:star6,
            markersize=8,
        )
    end

    if !isnothing(obs)
        _o = get_x_from_State(obs)
        x, y, z = split_vec_of_arrays(_o)
        scatter!(
            x, y, z,
            markershape=:circle,
            markerα = 0.5,
            #markersize=2,
            #markercolor="black"
        )
    end

    x_max = 1.0
    x_min = -1.0
    y_max = 0.2
    y_min = -1.0
    z_max = 2.0
    z_min = 0.0
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)*0.5
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2

    plot!(
        fig,
        xlims=(x_mid-max_range, x_mid+max_range),
        ylims=(y_mid-max_range, y_mid+max_range),
        zlims=(z_mid-max_range, z_mid+max_range),
        legend = false,
    )

    id = findall(x->x==true, jl)
    

    fig2 = plot(fig, camera=(90, 0))
    fig3 = plot(fig, camera=(0, 90))
    fig4 = plot(fig, camera=(0,0))
    fig_all = plot(
        fig, fig2, fig3, fig4, layout=(2, 2),
        size=(1000, 1000),
        title = string(round(t, digits=2)) * "[s]" * ", " * string(id)
    )
    return fig_all
end


"""アニメ制作"""
function make_animation(data, path)
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
        # _fig = draw_arm(
        #     data.q[i], data.dq[i], data.goal[i], data.obs[i], data.t[i], data.jl[i]
        # )

        _fig = draw_arm(
            data.nodes[i], data.goal[i], data.obs[i], data.t[i], data.jl[i]
        )

        frame(anim, _fig)
    end


    fname = path *  "animation.gif"

    gif(anim, fname, fps = 60)
    
    println("アニメ作成完了")
end 


"""ロボットアームを描画（nodesから）"""
function draw_arm(nodes, goal=nothing, obs=nothing, t=nothing, jl=nothing)


    fig = plot()
    x, y, z = nodes |> get_joint_positions_from_nodes |> split_vec_of_arrays
    plot!(
        fig,
        x, y, z,
        aspect_ratio = 1,
        marker=:circle,
        markerα = 0.5,
        label = "joints",
        xlabel = "X[m]", ylabel = "Y[m]", zlabel = "Z[m]"
    )

    cname = ("1", "2", "3", "4", "5", "6", "7", "GL")

    cpoints_x_global = nodes |> get_cpoints_from_nodes
    for (i, cs) in enumerate(cpoints_x_global)
        x, y, z = split_vec_of_arrays(cs)
        scatter!(
            fig,
            x, y, z,
            label = cname[i],
            markersize=5,
        )
    end

    if !isnothing(goal)
        scatter!(
            [goal.x[1]], [goal.x[2]], [goal.x[3]],
            markershape=:star6,
            markersize=8,
        )
    end

    if !isnothing(obs)
        _o = get_x_from_State(obs)
        x, y, z = split_vec_of_arrays(_o)
        scatter!(
            x, y, z,
            markershape=:circle,
            markerα = 0.5,
            #markersize=2,
            #markercolor="black"
        )
    end

    x_max = 1.0
    x_min = -1.0
    y_max = 0.2
    y_min = -1.0
    z_max = 2.0
    z_min = 0.0
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)*0.5
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2

    plot!(
        fig,
        xlims=(x_mid-max_range, x_mid+max_range),
        ylims=(y_mid-max_range, y_mid+max_range),
        zlims=(z_mid-max_range, z_mid+max_range),
        legend = false,
    )

    id = findall(x->x==true, jl)
    

    fig2 = plot(fig, camera=(90, 0))
    fig3 = plot(fig, camera=(0, 90))
    fig4 = plot(fig, camera=(0,0))
    fig_all = plot(
        fig, fig2, fig3, fig4, layout=(2, 2),
        size=(1000, 1000),
        title = string(round(t, digits=2)) * "[s]" * ", " * string(id)
    )
    return fig_all
end