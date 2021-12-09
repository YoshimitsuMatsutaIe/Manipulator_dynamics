
"""
siceで使用する平面冗長ロボットアーム
"""
module SiceDynamics

export calc_torque
export  calc_real_ddq

using LinearAlgebra

const m1 = 1.0
const m2 = 1.0
const m3 = 1.0
const m4 = 1.0

const l1 = 1.0
const l2 = 1.0
const l3 = 1.0
const l4 = 1.0

const lg1 = 0.5
const lg2 = 0.5
const lg3 = 0.5
const lg4 = 0.5

const g = 9.81

const I1 = 1.0
const I2 = 1.0
const I3 = 1.0
const I4 = 1.0


"""慣性行列"""
function M(q::Vector{T}) where T
    m_1_1 = I1 + I2 + I3 + I4 + l1^2*m2*sin(q[1])^2 + l1^2*m2*cos(q[1])^2 + l1^2*m3*sin(q[1])^2 + l1^2*m3*cos(q[1])^2 + l1^2*m4*sin(q[1])^2 + l1^2*m4*cos(q[1])^2 + 2.0*l1*l2*m3*sin(q[1])*sin(q[1] + q[2]) + 2.0*l1*l2*m3*cos(q[1])*cos(q[1] + q[2]) + 2.0*l1*l2*m4*sin(q[1])*sin(q[1] + q[2]) + 2.0*l1*l2*m4*cos(q[1])*cos(q[1] + q[2]) + 2.0*l1*l3*m4*sin(q[1])*sin(q[1] + q[2] + q[3]) + 2.0*l1*l3*m4*cos(q[1])*cos(q[1] + q[2] + q[3]) + 2.0*l1*lg2*m2*sin(q[1])*sin(q[1] + q[2]) + 2.0*l1*lg2*m2*cos(q[1])*cos(q[1] + q[2]) + 2.0*l1*lg3*m3*sin(q[1])*sin(q[1] + q[2] + q[3]) + 2.0*l1*lg3*m3*cos(q[1])*cos(q[1] + q[2] + q[3]) + 2.0*l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2^2*m3*sin(q[1] + q[2])^2 + l2^2*m3*cos(q[1] + q[2])^2 + l2^2*m4*sin(q[1] + q[2])^2 + l2^2*m4*cos(q[1] + q[2])^2 + 2.0*l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg1^2*m1*sin(q[1])^2 + lg1^2*m1*cos(q[1])^2 + lg2^2*m2*sin(q[1] + q[2])^2 + lg2^2*m2*cos(q[1] + q[2])^2 + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_1_2 = I2 + I3 + I4 + l1*l2*m3*sin(q[1])*sin(q[1] + q[2]) + l1*l2*m3*cos(q[1])*cos(q[1] + q[2]) + l1*l2*m4*sin(q[1])*sin(q[1] + q[2]) + l1*l2*m4*cos(q[1])*cos(q[1] + q[2]) + l1*l3*m4*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*l3*m4*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg2*m2*sin(q[1])*sin(q[1] + q[2]) + l1*lg2*m2*cos(q[1])*cos(q[1] + q[2]) + l1*lg3*m3*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*lg3*m3*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2^2*m3*sin(q[1] + q[2])^2 + l2^2*m3*cos(q[1] + q[2])^2 + l2^2*m4*sin(q[1] + q[2])^2 + l2^2*m4*cos(q[1] + q[2])^2 + 2.0*l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg2^2*m2*sin(q[1] + q[2])^2 + lg2^2*m2*cos(q[1] + q[2])^2 + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_1_3 = I3 + I4 + l1*l3*m4*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*l3*m4*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg3*m3*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*lg3*m3*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_1_4 = I4 + l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_2_1 = I2 + I3 + I4 + l1*l2*m3*sin(q[1])*sin(q[1] + q[2]) + l1*l2*m3*cos(q[1])*cos(q[1] + q[2]) + l1*l2*m4*sin(q[1])*sin(q[1] + q[2]) + l1*l2*m4*cos(q[1])*cos(q[1] + q[2]) + l1*l3*m4*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*l3*m4*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg2*m2*sin(q[1])*sin(q[1] + q[2]) + l1*lg2*m2*cos(q[1])*cos(q[1] + q[2]) + l1*lg3*m3*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*lg3*m3*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2^2*m3*sin(q[1] + q[2])^2 + l2^2*m3*cos(q[1] + q[2])^2 + l2^2*m4*sin(q[1] + q[2])^2 + l2^2*m4*cos(q[1] + q[2])^2 + 2.0*l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg2^2*m2*sin(q[1] + q[2])^2 + lg2^2*m2*cos(q[1] + q[2])^2 + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_2_2 = I2 + I3 + I4 + l2^2*m3*sin(q[1] + q[2])^2 + l2^2*m3*cos(q[1] + q[2])^2 + l2^2*m4*sin(q[1] + q[2])^2 + l2^2*m4*cos(q[1] + q[2])^2 + 2.0*l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + 2.0*l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + 2.0*l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg2^2*m2*sin(q[1] + q[2])^2 + lg2^2*m2*cos(q[1] + q[2])^2 + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_2_3 = I3 + I4 + l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_2_4 = I4 + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_3_1 = I3 + I4 + l1*l3*m4*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*l3*m4*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg3*m3*sin(q[1])*sin(q[1] + q[2] + q[3]) + l1*lg3*m3*cos(q[1])*cos(q[1] + q[2] + q[3]) + l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_3_2 = I3 + I4 + l2*l3*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*l3*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg3*m3*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3]) + l2*lg3*m3*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3]) + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_3_3 = I3 + I4 + l3^2*m4*sin(q[1] + q[2] + q[3])^2 + l3^2*m4*cos(q[1] + q[2] + q[3])^2 + 2.0*l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + 2.0*l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg3^2*m3*sin(q[1] + q[2] + q[3])^2 + lg3^2*m3*cos(q[1] + q[2] + q[3])^2 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_3_4 = I4 + l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_4_1 = I4 + l1*lg4*m4*sin(q[1])*sin(q[1] + q[2] + q[3] + q[4]) + l1*lg4*m4*cos(q[1])*cos(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_4_2 = I4 + l2*lg4*m4*sin(q[1] + q[2])*sin(q[1] + q[2] + q[3] + q[4]) + l2*lg4*m4*cos(q[1] + q[2])*cos(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_4_3 = I4 + l3*lg4*m4*sin(q[1] + q[2] + q[3])*sin(q[1] + q[2] + q[3] + q[4]) + l3*lg4*m4*cos(q[1] + q[2] + q[3])*cos(q[1] + q[2] + q[3] + q[4]) + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    m_4_4 = I4 + lg4^2*m4*sin(q[1] + q[2] + q[3] + q[4])^2 + lg4^2*m4*cos(q[1] + q[2] + q[3] + q[4])^2
    
    [
        m_1_1 m_1_2 m_1_3 m_1_4
        m_2_1 m_2_2 m_2_3 m_2_4
        m_3_1 m_3_2 m_3_3 m_3_4
        m_4_1 m_4_2 m_4_3 m_4_4
    ]
end


"""コリオリ項と重力項の和"""
function C_and_G(q::Vector{T}, dq::Vector{T}) where T
    u_all1 = g*lg1*m1*cos(q[1]) + g*m2*(l1*cos(q[1]) + lg2*cos(q[1] + q[2])) + g*m3*(l1*cos(q[1]) + l2*cos(q[1] + q[2]) + lg3*cos(q[1] + q[2] + q[3])) + g*m4*(l1*cos(q[1]) + l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) + lg4*cos(q[1] + q[2] + q[3] + q[4])) - 0.5*m2*((-2*dq[1]*l1*sin(q[1]) - 2*lg2*(dq[1] + dq[2])*sin(q[1] + q[2]))*(dq[1]*l1*cos(q[1]) + lg2*(dq[1] + dq[2])*cos(q[1] + q[2])) + (-dq[1]*l1*sin(q[1]) - lg2*(dq[1] + dq[2])*sin(q[1] + q[2]))*(-2*dq[1]*l1*cos(q[1]) - 2*lg2*(dq[1] + dq[2])*cos(q[1] + q[2]))) + 0.5*m2*((-2*l1*sin(q[1]) - 2*lg2*sin(q[1] + q[2]))*(-dq[1]^2*l1*cos(q[1]) - lg2*(dq[1] + dq[2])^2*cos(q[1] + q[2])) + (2*l1*cos(q[1]) + 2*lg2*cos(q[1] + q[2]))*(-dq[1]^2*l1*sin(q[1]) - lg2*(dq[1] + dq[2])^2*sin(q[1] + q[2])) + (-2*dq[1]*l1*sin(q[1]) - 2*lg2*(dq[1] + dq[2])*sin(q[1] + q[2]))*(dq[1]*l1*cos(q[1]) + lg2*(dq[1] + dq[2])*cos(q[1] + q[2])) + (-dq[1]*l1*sin(q[1]) - lg2*(dq[1] + dq[2])*sin(q[1] + q[2]))*(-2*dq[1]*l1*cos(q[1]) - 2*lg2*(dq[1] + dq[2])*cos(q[1] + q[2]))) - 0.5*m3*((-2*dq[1]*l1*sin(q[1]) - 2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3])) + (-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*(-2*dq[1]*l1*cos(q[1]) - 2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]))) + 0.5*m3*((-2*l1*sin(q[1]) - 2*l2*sin(q[1] + q[2]) - 2*lg3*sin(q[1] + q[2] + q[3]))*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3])) + (2*l1*cos(q[1]) + 2*l2*cos(q[1] + q[2]) + 2*lg3*cos(q[1] + q[2] + q[3]))*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3])) + (-2*dq[1]*l1*sin(q[1]) - 2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3])) + (-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*(-2*dq[1]*l1*cos(q[1]) - 2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]))) - 0.5*m4*((-2*dq[1]*l1*sin(q[1]) - 2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])) + (-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(-2*dq[1]*l1*cos(q[1]) - 2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))) + 0.5*m4*((-2*l1*sin(q[1]) - 2*l2*sin(q[1] + q[2]) - 2*l3*sin(q[1] + q[2] + q[3]) - 2*lg4*sin(q[1] + q[2] + q[3] + q[4]))*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*cos(q[1] + q[2] + q[3] + q[4])) + (2*l1*cos(q[1]) + 2*l2*cos(q[1] + q[2]) + 2*l3*cos(q[1] + q[2] + q[3]) + 2*lg4*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*sin(q[1] + q[2] + q[3] + q[4])) + (-2*dq[1]*l1*sin(q[1]) - 2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])) + (-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(-2*dq[1]*l1*cos(q[1]) - 2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])))
    u_all2 = g*lg2*m2*cos(q[1] + q[2]) + g*m3*(l2*cos(q[1] + q[2]) + lg3*cos(q[1] + q[2] + q[3])) + g*m4*(l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) + lg4*cos(q[1] + q[2] + q[3] + q[4])) - 0.5*m2*(-2*lg2*(dq[1] + dq[2])*(-dq[1]*l1*sin(q[1]) - lg2*(dq[1] + dq[2])*sin(q[1] + q[2]))*cos(q[1] + q[2]) - 2*lg2*(dq[1] + dq[2])*(dq[1]*l1*cos(q[1]) + lg2*(dq[1] + dq[2])*cos(q[1] + q[2]))*sin(q[1] + q[2])) + 0.5*m2*(-2*lg2*(dq[1] + dq[2])*(-dq[1]*l1*sin(q[1]) - lg2*(dq[1] + dq[2])*sin(q[1] + q[2]))*cos(q[1] + q[2]) - 2*lg2*(dq[1] + dq[2])*(dq[1]*l1*cos(q[1]) + lg2*(dq[1] + dq[2])*cos(q[1] + q[2]))*sin(q[1] + q[2]) + 2*lg2*(-dq[1]^2*l1*sin(q[1]) - lg2*(dq[1] + dq[2])^2*sin(q[1] + q[2]))*cos(q[1] + q[2]) - 2*lg2*(-dq[1]^2*l1*cos(q[1]) - lg2*(dq[1] + dq[2])^2*cos(q[1] + q[2]))*sin(q[1] + q[2])) - 0.5*m3*((-2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3])) + (-2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]))*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))) + 0.5*m3*((-2*l2*sin(q[1] + q[2]) - 2*lg3*sin(q[1] + q[2] + q[3]))*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3])) + (2*l2*cos(q[1] + q[2]) + 2*lg3*cos(q[1] + q[2] + q[3]))*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3])) + (-2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3])) + (-2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]))*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))) - 0.5*m4*((-2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])) + (-2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))) + 0.5*m4*((-2*l2*sin(q[1] + q[2]) - 2*l3*sin(q[1] + q[2] + q[3]) - 2*lg4*sin(q[1] + q[2] + q[3] + q[4]))*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*cos(q[1] + q[2] + q[3] + q[4])) + (2*l2*cos(q[1] + q[2]) + 2*l3*cos(q[1] + q[2] + q[3]) + 2*lg4*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*sin(q[1] + q[2] + q[3] + q[4])) + (-2*l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])) + (-2*l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - 2*l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4])))
    u_all3 = g*lg3*m3*cos(q[1] + q[2] + q[3]) + g*m4*(l3*cos(q[1] + q[2] + q[3]) + lg4*cos(q[1] + q[2] + q[3] + q[4])) - 0.5*m3*(-2*lg3*(dq[1] + dq[2] + dq[3])*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*cos(q[1] + q[2] + q[3]) - 2*lg3*(dq[1] + dq[2] + dq[3])*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]))*sin(q[1] + q[2] + q[3])) + 0.5*m3*(-2*lg3*(dq[1] + dq[2] + dq[3])*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]))*cos(q[1] + q[2] + q[3]) - 2*lg3*(dq[1] + dq[2] + dq[3])*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + lg3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]))*sin(q[1] + q[2] + q[3]) + 2*lg3*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3]))*cos(q[1] + q[2] + q[3]) - 2*lg3*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - lg3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3]))*sin(q[1] + q[2] + q[3])) - 0.5*m4*((-2*l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])) + (-2*l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))) + 0.5*m4*((-2*l3*sin(q[1] + q[2] + q[3]) - 2*lg4*sin(q[1] + q[2] + q[3] + q[4]))*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*cos(q[1] + q[2] + q[3] + q[4])) + (2*l3*cos(q[1] + q[2] + q[3]) + 2*lg4*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*sin(q[1] + q[2] + q[3] + q[4])) + (-2*l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])) + (-2*l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4])))
    u_all4 = g*lg4*m4*cos(q[1] + q[2] + q[3] + q[4]) - 0.5*m4*(-2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*cos(q[1] + q[2] + q[3] + q[4]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))*sin(q[1] + q[2] + q[3] + q[4])) + 0.5*m4*(-2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*(-dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]))*cos(q[1] + q[2] + q[3] + q[4]) - 2*lg4*(dq[1] + dq[2] + dq[3] + dq[4])*(dq[1]*l1*cos(q[1]) + l2*(dq[1] + dq[2])*cos(q[1] + q[2]) + l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) + lg4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]))*sin(q[1] + q[2] + q[3] + q[4]) + 2*lg4*(-dq[1]^2*l1*sin(q[1]) - l2*(dq[1] + dq[2])^2*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*sin(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*sin(q[1] + q[2] + q[3] + q[4]))*cos(q[1] + q[2] + q[3] + q[4]) - 2*lg4*(-dq[1]^2*l1*cos(q[1]) - l2*(dq[1] + dq[2])^2*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])^2*cos(q[1] + q[2] + q[3]) - lg4*(dq[1] + dq[2] + dq[3] + dq[4])^2*cos(q[1] + q[2] + q[3] + q[4]))*sin(q[1] + q[2] + q[3] + q[4]))

    [
        u_all1
        u_all2
        u_all3
        u_all4
    ]
end



"""トルクを計算

計算トルク法です  
q : 関節角度ベクトル  
dq : 関節角速度ベクトル  
desired_ddq : （所望の）関節角速度ベクトル  
"""
function calc_torque(
    q::Vector{TU}, dq::Vector{TU}, desired_ddq::Vector{TU}
    ) where TU
    M(q)*desired_ddq .+ C_and_G(q, dq) |> vec
end


# """物体から受ける力を計算"""
# function calc_Fc(;
#     u::Vector{TU}, q::Vector{TU}, dq::Vector{TU},
#     F::Vector{TU},
#     Jend::Matrix{TU}
#     ) where TU

#     inv(transpose(Jend)) * (M(q)*desired_ddq .+ C_and_G(q, dq) .- u .- F) |> vec
# end

"""現実世界での加速度

u : 入力トルクベクトル R(7)  
q : （現在の）関節角度ベクトル  
dq : （現在の）関節角速度ベクトル  
F : 外乱ベクトル R(7)  
Fc : エンドエフェクタに加わる外力ベクトル  
Jend : エンドエフェクタのヤコビ行列  
"""
function calc_real_ddq(;
    u::Vector{TU}, q::Vector{TU}, dq::Vector{TU},
    F::Vector{TU}, Fc::Vector{TU},
    Jend::Matrix{TU}
    ) where TU

    # println(size(transpose(Jend)))
    # println(transpose(Jend) * Fc)
    # println(F)

    _Fc = transpose(Jend) * Fc
    @. F += _Fc

    # println(F)

    inv(M(q)) * (u .+ F .- (C_and_G(q, dq))) |> vec
end

end

# function test()
#     q = [
#         -90.0
#         0.0
#         0.0
#         0.0
#     ] * pi / 180
#     dq = zero(q)
#     ddq = zero(q)


#     println(M(q))

#     @time println(calc_torque(q, dq, ddq))
# end


# test()