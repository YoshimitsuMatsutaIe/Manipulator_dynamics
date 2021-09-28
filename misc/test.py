
r_bar_1 = [
    np.array([[0, L1/2, -L0/2, 1]]).T,
    np.array([[0, -L1/2, -L0/2, 1]]).T,
    np.array([[L1/2, 0, -L0/2, 1]]).T,
    np.array([[-L1/2, 0, -L0/2, 1]]).T,
]  # 1座標系からみた制御点位置

r_bar_2 = [
    np.array([[0, 0, L3/2, 1]]).T,
    np.array([[0, 0, -L3/2, 1]]).T,
]

r_bar_3 = [
    np.array([[0, L3/2, -L2*2/3, 1]]).T,
    np.array([[0, -L3/2, -L2*2/3, 1]]).T,
    np.array([[L3/2, 0, -L2*2/3, 1]]).T,
    np.array([[-L3/2, 0, -L2*2/3, 1]]).T,
    np.array([[0, L3/2, -L2*1/3, 1]]).T,
    np.array([[0, -L3/2, -L2*1/3, 1]]).T,
    np.array([[L3/2, 0, -L2*1/3, 1]]).T,
    np.array([[-L3/2, 0, -L2*1/3, 1]]).T,
]

r_bar_4 = [
    np.array([[0, 0, L3/2, 1]]).T,
    np.array([[0, 0, -L3/2, 1]]).T,
]

r_bar_5 = [
    np.array([[0, L5/2, -L4/2, 1]]).T,
    np.array([[0, -L5/2, -L4/2, 1]]).T,
    np.array([[L5/2, 0, -L4/2, 1]]).T,
    np.array([[-L5/2, 0, -L4/2, 1]]).T,
]

r_bar_6 = [
    np.array([[0, 0, L5/2, 1]]).T,
    np.array([[0, 0, -L5/2, 1]]).T,
]

r_bar_7 = [
    np.array([[0, L5/2, L6/2, 1]]).T,
    np.array([[0, -L5/2, L6/2, 1]]).T,
    np.array([[L5/2, 0, L6/2, 1]]).T,
    np.array([[-L5/2, 0, L6/2, 1]]).T,
]

r_bar_GL = [
    np.array([[0, 0, 0, 1]]).T
]

# 追加
r_bars = [
    r_bar_1, r_bar_2, r_bar_3, r_bar_4, r_bar_5, r_bar_6, r_bar_7, r_bar_GL,
]
