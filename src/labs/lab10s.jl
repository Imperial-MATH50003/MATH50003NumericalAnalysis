# # MATH50003 (2024–25)
# # Lab 10: VI.3 Gaussian Quadrature

# **Learning Outcomes**
#
# Mathematical knowledge:


# Coding knowledge:

# 1. `eigen` for computing eigenvalues and eigenvectors

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


# **Problem 1(a)** Compute all 100 roots of the Chebyshev U polynomial $U_{100}(x)$.

## TODO: Use a truncation Jacobi matrix to compute the roots of U₁₀₀(x)

# **Problem 1(b)** Compute all 100 roots of the Legendre polynomial $P_{100}(x)$.

## TODO: Use a truncation Jacobi matrix to compute the roots of P₁₀₀(x)

# ## VI.3.2 Properties of Gaussian quadrature