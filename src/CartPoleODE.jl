module CartPoleODE

using StaticArrays

const nx = 4    # number of states
const nu = 1    # number of inputs

struct Model
    g     # ms⁻² - gravity
    m_c   # kg   - mass of the cart
    m_p   # kg   - mass of the point-mass 
    l     # m    - length of the pole
    b_c   # Ns⁻¹ - damping coefficient for the cart
    b_p   # Ns⁻¹ - damping coefficient for the pole
end

mass_matrix(m, x) = @SMatrix [
    m.m_c+m.m_p m.m_p*m.l*cos(x[2])
    m.m_p*m.l*cos(x[2]) m.m_p*m.l^2
]

torque_vector(m, x, u) = @SVector [
    m.m_p * m.l * sin(x[2]) * x[4]^2 + u[1] - m.b_c * x[2],
    -m.g * m.m_p * m.l * sin(x[2]) - m.b_p * x[4]
]

function f(m, x, u)
    M = mass_matrix(m, x)
    τ = torque_vector(m, x, u)

    return vcat(x[3:4], M \ τ)
end

function f!(m, ẋ, x, u)
    M = mass_matrix(m, x)
    τ = torque_vector(m, x, u)

    @views ẋ[1:2] .= x[3:4]
    @views ẋ[3:4] .= M \ τ
end

end # module CartPoleODE
