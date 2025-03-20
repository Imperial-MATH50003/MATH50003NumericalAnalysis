# # MATH50003 (2024–25)
# # Lab 10: VI.3 Gaussian Quadrature

# We complete the module where we started: computing integrals via quadrature. 
# We saw how interpolation leads to interpolatory quadrature rules, but we were
# left with the question on how to choose the points. Orthogonal polynomials
# provide the answer: The roots (zeros) of degree $n$ polynomials give interpolatory
# quadrature rules with the special property that they are exact for polynomials of twice
# the degree. Moreover, they lead to effective transforms from values of functions to
# coefficients in an orthogonal polynomial expansion. In this lab we realise these
# results numerically. 

# **Learning Outcomes**
#
# Mathematical knowledge:

# 1. Roots of orthogonal polynomials and truncated Jacobi matrices.
# 2. Gaussian quadrature via eigenvalues of a truncated Jacobi matrix.
# 3. Transforms via Gaussian quadrature.

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

# We can see this example with Chebyshev polynomials. We first construct a $10×10$ truncation of the Jacobi matrix
# and compute its eigenvalues:

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


# **Problem 1(b)** Compute all 100 roots of the Legendre polynomial $P_{100}(x)$
# by constructing the multiplication matrix, symmetrising it and computing its eigenvalues.

## TODO: Use a truncation Jacobi matrix to compute the roots of P₁₀₀(x).


# ------

# ## VI.3.2 Properties of Gaussian quadrature

# Gaussian quadrature is a special quadrature rule that is exact for almost twice the degree of polynomials
# as other interpolatory quadrature rules. It is defined in terms of the eigenvalues and eigenvectors
# of the truncated Jacobi matrix as follows, using the Chebyshev T polynomials:

function gaussquadrature(J::SymTridiagonal)
    x,Q = eigen(J) # eigen computes the eigenvalues and eigenvectors.
    x, Q[1,:].^2
end


function gausschebyshevt(n)
    J = SymTridiagonal(zeros(n), [1/sqrt(2); fill(1/2, n-2)]) # symmetric Jacobi matrix
    x,w = gaussquadrature(J)
    x, π*w # ∫dx/sqrt(1-x^2)= π
end

x,w = gausschebyshevt(15)
@test w'exp.(x) ≈ 3.977463260506422 # = ∫exp(x)*dx/sqrt(1-x^2)

# This approximation actually converges faster than exponentially!

plot(2:15, [((x,w) = gausschebyshevt(n); abs(w'exp.(x) - 3.977463260506422)) for n=2:15]; yscale=:log10, label="error", title="Gauss–Chebyshev Quadrature")

# Proving the convergence for general Gaussian quadrature requires complex analysis and potential theory results that are a bit too
# advanced for Year 2. But in the specific case of Gauss–Chebyshev quadrature it can be proven via the connection with
# the Cosine expansion.


# We can also  use Gauss quadrature to compute expansion coefficients. To do this
# we need to construct the orthonormal polynomials $q_0(x) = 1/\sqrt{π}$ and $q_n(x) = T_n(x)\sqrt{2/π}$:

function orthonormalchebyshevt(n, x) # normalized T_n(x), 
    if n == 0
        1/sqrt(π)
    else
        cos(n*acos(x))*sqrt(2/π)
    end
end

# The inverse transform from coefficients to values is the Vandermonde-like matrix:
# $$
# V_n = \begin{bmatrix}
# q_0(x_1) & ⋯ & q_{n-1}(x_1) \\
# ⋮ & ⋱ & ⋱ \\
# q_0(x_n) & ⋯ & q_{n-1}(x_n)
# \end{bmatrix}
# $$
# The transform from values to coefficients of $q_n(x)$ is then $V_n^{-1}$, but we have a simple formula:
# $$
# V_n^{-1} = V_n^⊤{\rm diag}(w_1,…,w_n).
# $$
# We implement this as follows:

function chebyshevttransform(n)
    x,w = gausschebyshevt(n)
    V = orthonormalchebyshevt.((0:n-1)', x) # Vandermonde-like matrix for normalized T_n(x)
    V'*Diagonal(w) # transform matrix
end

n = 15
Q = chebyshevttransform(n)
c = Q*exp.(x) # expansion coefficients
@test c'*orthonormalchebyshevt.((0:n-1), 0.1) ≈ exp(0.1)

# We can see that the above actually interpolates, using a more
# oscillatory example where it's visually distinctive:

f = x -> cos(20x)
c = Q*f.(x) # expansion coefficients
g = range(-1,1,1000)
plot(g, f.(g); label="f")
plot!(g, orthonormalchebyshevt.((0:n-1)', g)*c; label="fₙ")
scatter!(x, f.(x); label=nothing)

# Taking a few more points and the approximation will have converged to high accuracy,
# actually super-exponentially fast.

# ---
# **Problem 2** Compute the Gauss–Chebyshev U quadrature rule for $w(x) = \sqrt{1-x^2}$ on $[-1,1]$.

function gausschebyshevu(n)
    ## TODO: implement Gauss–Chebyshev U quadrature
    
end

x,w = gausschebyshevu(3)
@test x ≈ [-1/sqrt(2),0,1/sqrt(2)]
@test w ≈ [0.3926990816987255, 0.7853981633974466, 0.3926990816987243]


# **Problem 3(a)** Compute the Gauss–Legendre quadrature rule, for $w(x) = 1$ on $[-1,1]$.

function gausslegendre(n)
    ## TODO: Compute the Gauss–Legendre quadrature rule for a uniform weight.
    
end
x,w = gausslegendre(3)
@test x ≈ [-sqrt(3/5), 0, sqrt(3/5)]
@test w ≈ [5/9,8/9,5/9]

n = 5
x,w = gausslegendre(n)
@test w'exp.(x) ≈ exp(1)-exp(-1) # even just 5 points converges

## We saw in Theorem 20 that it is exact for about twice the degree
## of polynomial. We can test this numerically:

for k = 0:2n-1
    if iseven(k)
        @test w'x.^k ≈ 2/(k+1)
    else
        @test abs(w'x.^k) ≤ 1E-14
    end
end

# **Problem 3(b)** Implement the following function that computes `orthonormallegendrep(n,x)` corresponding
# to the orthonormalised Legendre polynomials.

function orthonormallegendrep(n, x)
    ## TODO: implement the orthonormalised Legendre polynomials
    
end

@test orthonormallegendrep(5, 0.1) ≈ 0.41939059365476206

# **Problem 3(c)** Implement the Legendre transform from values of a function on Gauss–Legendre points to coefficients in
# a Legendre expansion.

function legendreptransform(n)
    ## TODO: Construct the n × n matrix mapping from samples at zeros of Legendre polynonials to coefficients
    
end

n = 15
x,w = gausslegendre(n)
Q = legendreptransform(n)
c = Q*exp.(x) # expansion coefficients
@test c'*orthonormallegendrep.((0:n-1), 0.1) ≈ exp(0.1)
