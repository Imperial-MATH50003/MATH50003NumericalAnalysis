# # MATH50003 (2024–25)
# # Lab 10: VI.3 Gaussian Quadrature

# **Learning Outcomes**
#
# Mathematical knowledge:


# Coding knowledge:

# 1. `eigvals` and `eigen` for computing eigenvalues and eigenvectors

using Plots, LinearAlgebra

# ## VI.3.1 Roots of orthogonal polynomials and truncated Jacobi matrices

# Gaussian quadrature is defined in terms of the roots of orthogonal polynomials.
# We know how to compute orthogonal polynomials so in theory we could use Newton's method
# (with dual numbers to compute the derivatives), however, we need _all_ roots and
# it is a challenging problem choosing initial guesses.

# Instead, we will recast the roots as eigenvalues of a symmetric tridiagonal matrix.
# Having previously thought of eigenvalues as roots of a characteristic polynomial this
# may appear to be a circular idea. But it turns out that there are very effective
# algorithms for computing all eigenvalues (which are taught in Computational Linear Algebra),
# and in fact finding roots of polynomials is typically done via eigenvalue algorithms.
# For us we can use the `eigvals` function to compute these.

# We can see this example with Chebyshev polynomials. We first construct a $10×10$ truncation of the Jacobi matrix:

n = 10
J = SymTridiagonal(zeros(n), [1/sqrt(2); fill(1/2, n-2)])
x = eigvals(J)

# By plotting we can see that we have successfully computed
# the roots of $T_{10}(x)$:

g = range(-1,1,100)
plot(g, cos.(n*acos.(g)); label="T₁₀(x)")
scatter!(x, zero(x); label="roots")

# -----

# **Problem 1(a)** Compute all 50 roots of the Chebyshev U polynomial $U_{50}(x)$.

## TODO: Use a truncation Jacobi matrix to compute the roots of U₅₀(x)
## SOLUTION
## The Chebyshev U polynomials are already orthonormal with a simple Jacobi matrix:
g = range(-1,1,1000)
n = 50
J = SymTridiagonal(zeros(n), fill(1/2, n-1))
x = eigvals(J)

plot(g, sin.((n+1)*acos.(g)) ./ sin.(acos.(g)); label="U₁₀₀(x)", ylim=(-2,2))
scatter!(x, zero(x); label="roots")
## END

# **Problem 1(b)** Compute all 100 roots of the Legendre polynomial $P_{100}(x)$
# by constructing the multiplication matrix, symmetrising it and computing its eigenvalues.

## TODO: Use a truncation Jacobi matrix to compute the roots of P₁₀₀(x).
## SOLUTION
## This is harder since the Jacobi matrix is not symmetric. But we can
## adapt the symmetrise code from last lecture to handle non-monic polynomials:
function symmetrise(X::Tridiagonal)
    n = size(X,1)
    k = zeros(n)
    k[1] = 1 # The normalisation for q_0(x) doesn't impact the Jacobi matrix
    for j = 1:n-1
        k[j+1] = k[j]*sqrt(X[j+1,j]/X[j,j+1])
    end
    SymTridiagonal(X.d,  X.dl .* k[1:n-1] ./ k[2:n])
end

n = 100
X = Tridiagonal((1:(n-1)) ./ (1:2:(2n-3)), zeros(n), (1:(n-1)) ./ (3:2:(2n-1)))
x = eigvals(symmetrise(X))

function legendrep(n, x)
    Pₖ₋₁ = 1.0
    if n == 0
        return Pₖ₋₁
    end
    Pₖ = x
    for k = 1:n-1
        Pₖ,Pₖ₋₁ = (2k+1)/(k+1) * x*Pₖ - k/(k+1)*Pₖ₋₁, Pₖ
    end
    Pₖ
end

g = range(-1,1,2000)
plot(g,legendrep.(n, g); ylims=(-0.3,0.3), label="P₁₀₀")
scatter!(x, zero(x); label="roots")

## END

# ------

# ## VI.3.2 Properties of Gaussian quadrature

# Gaussian quadrature is a special quadrature rule that is exact for twice the number of polynomials
# as other interpolatory quadrature rules. It is defined in terms of the eigenvalues and eigenvectors
# of the truncated Jacobi matrix as follows, using the Chebyshev T polynomials:

function gaussquadrature(J::SymTridiagonal)
    x,Q = eigen(J) # eigen computes the eigenvalues and eigenvectors.
    x, Q[1,:].^2
end

n = 15
J = SymTridiagonal(zeros(n), [1/sqrt(2); fill(1/2, n-2)])
x,w = gaussquadrature(J)
μ = π # = ∫ dx/sqrt(1-x^2)
@test μ*w'exp.(x) ≈ 3.977463260506422 # = ∫exp(x)*dx/sqrt(1-x^2)

# We can use gauss quadrature to compute expansion coefficients via:

function orthonormalchebyshevt(n, x) # normalized T_n(x), 
    if n == 0
        1/sqrt(π)
    else
        cos(n*acos(x))*sqrt(2/π)
    end
end

V = orthonormalchebyshevt.((0:n-1)', x) # Vandermonde-like matrix for normalized T_n(x)
Q = V'*Diagonal(μ*w) # transform matrix
c = Q*exp.(x) # expansion coefficients
@test c'*orthonormalchebyshevt.((0:n-1), 0.1) ≈ exp(0.1)

# We can see that the above actually interpolates, using an example
# where we can see it's not exact:

f = x -> cos(20x)
c = Q*f.(x) # expansion coefficients
g = range(-1,1,1000)
plot(g, f.(g); label="f")
plot!(g, orthonormalchebyshevt.((0:n-1)', g)*c; label="fₙ")
scatter!(x, f.(x); label=nothing)

# Taking a few more points and the approximation will have converged to high accuracy.

# ---
# **Problem 2(a)**