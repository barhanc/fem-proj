#=
Numerical solution of a certain type of ODE using Finite Element Method written
as a project for Differential Equations course at AGH UST.

Author: Bartosz Hanc
Git repository: https://github.com/Bh4n5/finite_element_method.git
=#

using Plots

# Numerical integration using 3-point Gauss-Legendre quadrature for every
# interval
function gaussQuad(f, a, b)
    c = (b - a) / 2
    g = x -> (b - a) / 2 * x + (a + b) / 2
    x = [-sqrt(3 / 5), 0, sqrt(3 / 5)]
    w = [5 / 9, 8 / 9, 5 / 9]
    return c * sum([w[i] * f(g(x[i])) for i = 1:3])
end

# Integral operator
function integral(f, a, b)
    n = 10^3
    dx = (b - a) / n
    return sum([gaussQuad(f, a + i * dx, a + (i + 1) * dx) for i = 0:n-1])
end

#= Main solver for eq. (E(x)u'(x))' = 0 with boundary conditions u(2) = 0, u'(0)
+ u(0) = 10 using discretization of weak formulation: 
∀i: ∑ aⱼ[ϕ(vᵢ,vⱼ)-E(0)vᵢ(0)vⱼ(0)]=-10E(0)vᵢ(0) ⟹ B⋅A=L where
B=[ϕ(vᵢ,vⱼ)-E(0)vᵢ(0)vⱼ(0)]ₙₓₙ , A=[aᵢ]ₙ and L=[-10E(0)vᵢ(0)]ₙ =#
function solve(n)
    # Input parameters
    a, b = 0, 2
    dx = (b - a) / n
    E = x -> x >= 0 && x <= 1 ? 3 : 5

    isIn = (k, x) -> x <= a + k * dx && x >= a + (k - 1) * dx
    y₊ = (k, x) -> 1 + (x - a - k * dx) / dx
    y₋ = (k, x) -> 1 - (x - a - k * dx) / dx

    #= Define families of hat functions (composed using linear functions defined
    above) and its derivatives (denoted with prefix 'D') used as a basis of a
    discrete space of test functions. =#
    v = (k, x) -> isIn(k, x) ? y₊(k, x) : isIn(k + 1, x) ? y₋(k, x) : 0
    Dv = (k, x) -> isIn(k, x) ? 1 / dx : isIn(k + 1, x) ? -1 / dx : 0

    #= Define map ϕ:F×F⟶ℜ which takes a family of functions f = {f₁,f₂, ...},
    where fᵢ ∈ F, indices i and j and returns ∫ₐᵇ E(x)fᵢfⱼ dx .=#
    ϕ = (f, i, j) -> integral(x -> E(x) * f(i, x) * f(j, x), a, b)

    # Define delta function δ:N×N⟶{0,1} as 1 if i=j or i-j=±1 | 0 otherwise
    δ = (i, j) -> i == j || abs(i - j) == 1

    # Define matrix B and vector L. 
    B::Matrix{Float64} = [δ(i, j) ? ϕ(Dv, i, j) - E(0) * v(i, 0) * v(j, 0) : 0 for i in 0:n-1, j in 0:n-1]
    L::Vector{Float64} = [-10 * E(0) * v(i, 0) for i in 0:n-1]

    # Solve the above system of equations for coefficients aᵢ using Julia's `\`
    # operator efficient for sparse matrices.
    A = B \ L

    # Solution
    return x -> sum([A[i+1] * v(i, x) for i in 0:n-1])
end

function IO()
    print(">> Input number of elements: ")
    n = parse(Int64, readline())

    u, x = solve(n), range(0, 2, 100)
    plt = plot(x, u.(x),
        title="Solution",
        label="u(x)",
        size=(600, 400)
    )
    display(plt)

    println(">> Press enter to exit...")
    readline()

    return nothing
end

IO()