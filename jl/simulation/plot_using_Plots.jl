"""データ可視化"""



using Plots


"""基本的なものだけplot"""
function plot_simulation_data(data)
    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.q)
    fig_q = plot(data.t, q1, ylabel="q", label="1")
    plot!(fig_q, data.t, q2, label="2")
    plot!(fig_q, data.t, q3, label="3")
    plot!(fig_q, data.t, q4, label="4")
    plot!(fig_q, data.t, q5, label="5")
    plot!(fig_q, data.t, q6, label="6")
    plot!(fig_q, data.t, q7, label="7")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.dq)
    fig_dq = plot(data.t, q1, ylabel="dq", label="1")
    plot!(fig_dq, data.t, q2, label="2")
    plot!(fig_dq, data.t, q3, label="3")
    plot!(fig_dq, data.t, q4, label="4")
    plot!(fig_dq, data.t, q5, label="5")
    plot!(fig_dq, data.t, q6, label="6")
    plot!(fig_dq, data.t, q7, label="7")

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.ddq)
    fig_ddq = plot(data.t, q1, ylabel="ddq", label="1")
    plot!(fig_ddq, data.t, q2, label="2")
    plot!(fig_ddq, data.t, q3, label="3")
    plot!(fig_ddq, data.t, q4, label="4")
    plot!(fig_ddq, data.t, q5, label="5")
    plot!(fig_ddq, data.t, q6, label="6")
    plot!(fig_ddq, data.t, q7, label="7")

    fig_error = plot(data.t, data.error, label="error", ylabel="error [m]", ylims=(0.0,))
    fig_dis_to_obs = plot(data.t, data.dis_to_obs, label="min dis to obs", ylabel="min distance to obs [m]", ylims=(0.0,))

    fig = plot(
        fig_q, fig_dq, fig_ddq, fig_error, fig_dis_to_obs,
        layout=(5,1),
        size=(500,1200)
    )

    return fig
end




function draw_arm(q=q_neutral, dq=zeros(Float64, 7), goal=nothing, obs=nothing, t=nothing)

    _, _, _, _, _, _, _, _, cpoints_x_global, _, joints_x_global, _, = calc_all(q, dq)
    #fig = plot(size=(800,700))

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


    fig2 = plot(fig, camera=(90, 0))
    fig3 = plot(fig, camera=(0, 90))
    fig4 = plot(fig, camera=(0,0))
    fig_all = plot(
        fig, fig2, fig3, fig4, layout=(2, 2),
        size=(1000, 1000),
        title = string(round(t, digits=2)) * "[s]"
    )
    return fig_all
end


"""アニメ制作"""
function make_animation(data)
    # 枚数決める
    epoch_max = 100
    epoch = length(data.t)
    if epoch < epoch_max
        step = 1
    else
        step = div(epoch, epoch_max)
    end

    println(step)

    anim = Animation()
    @gif for i in 1:step:length(data.q)
        _fig = draw_arm(
            data.q[i], data.dq[i], data.goal[i], data.obs[i], data.t[i]
        )
        frame(anim, _fig)
    end
    gif(anim, "test5.gif", fps = 60)
    #return anim
end 
