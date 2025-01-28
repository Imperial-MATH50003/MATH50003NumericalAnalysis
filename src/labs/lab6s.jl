# # MATH50003 (2024–25)
# # Lab 6: III.4 QR Factorisation and IV.1 Polynomial Regression


# We also compute the QR factorisation with Householder reflections, and use this
# to solve least squares problems.


# We also explore polynomial interpolation and regression, and see that when
# interpolating at an evenly spaced grid one can encounter issues with convergence.
# This is overcome via regression, but we are left with the question of how to
# solve the underlying least squares problems. 


# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 2. Computing the QR factorisation using reflections.
# 3. Computing a tridiagonal QR factorisation using rotations.
# 4. The relationship between QR and least squares.
# 2. Vandermonde matrices and least squares.
# 3. Issues with interpolation at evenly spaced points with functions with small radii of convergence.

#
# Coding knowledge:
#
# 1. The `atan(y,x)` function and the `I` convenience syntax.
# 2. Templating fields in a type.
# 2. Solving least squares problems via `\`.
# 3. Using the `qr` function to solve least squares.



using LinearAlgebra, Test





# -----


# ## III.6 QR Factorisation

# The QR factorisation of a matrix $A ∈ ℂ^{m × n}$ is of the form
# $$
# A = QR
# $$
# where $Q$ is unitary and $R$ is right-triangular: all entries below the diagonal are zero. We focus on the case where $m ≥ n$. 
# It can be computed using Gram–Schmidt, Householder reflections or rotations.

# ### III.6.1 Reduced QR and Gram–Schmidt

# The Gram–Schmidt process can be used to compute the QR factorisation by orthogonalising the columns
# of $A$ in sequence. We won't discuss this in more detail as it is numerically better to use reflections/rotations.

# ### III.6.2 Householder reflections and QR

# In the notes we use Householder reflections to prove that a QR factorisation exists.
# The iterative proof actually encodes an algorithm, which we can implement as follows:


function householderqr(A)
    T = eltype(A)
    m,n = size(A)
    if n > m
        error("More columns than rows is not supported")
    end

    R = zeros(T, m, n)
    Q = Matrix(one(T)*I, m, m)
    Aⱼ = copy(A)

    for j = 1:n
        𝐚₁ = Aⱼ[:,1] # first columns of Aⱼ
        Q₁ = dense_householderreflection(𝐚₁)
        Q₁Aⱼ = Q₁*Aⱼ # multiply Aⱼ by the Householder reflection
        α,𝐰 = Q₁Aⱼ[1,1],Q₁Aⱼ[1,2:end]

        ## populate returned data
        R[j,j] = α
        R[j,j+1:end] = 𝐰

        ## following is equivalent to Q = Q*[I 0 ; 0 Qⱼ]
        Q[:,j:end] = Q[:,j:end]*Q₁

        Aⱼ = Q₁Aⱼ[2:end,2:end] # this is the "induction"
    end
    Q,R
end

m,n = 100,50
A = randn(m,n)
Q,R = householderqr(A)
@test Q'Q ≈ I
@test Q*R ≈ A


# Note because we are forming a full matrix representation of each Householder
# reflection this is a slow algorithm: it uses $O(m^2 n^2)$ operations, which is too many!
# By being more careful about how we apply and store reflections we can avoid this,
# in particular, taking advantage of the types `Reflection` and `Reflections`.

# ------

# **Problem 3** Complete the following function that implements
# Householder QR for a real matrix $A ∈ ℝ^{m × n}$ where $m ≥ n$ using only $O(mn^2)$ operations, using 
#  `Reflection` and `Reflections`.

function householderqr(A)
    T = eltype(A)
    m,n = size(A)
    if n > m
        error("More columns than rows is not supported")
    end

    R = zeros(T, m, n)
    Q = Reflections(zeros(T, m, n))
    Aⱼ = copy(A)

    for j = 1:n
        ## TODO: rewrite householder QR to use Reflection,
        ## Reflections and householderreflection, in a way that one achieves O(mn^2) operations
        ## SOLUTION
        𝐚₁ = Aⱼ[:,1] # first columns of Aⱼ
        Q₁ = householderreflection(𝐚₁[1] < 0, 𝐚₁)
        Q₁Aⱼ = Q₁*Aⱼ
        α,𝐰 = Q₁Aⱼ[1,1],Q₁Aⱼ[1,2:end]
        Aⱼ₊₁ = Q₁Aⱼ[2:end,2:end]

        ## populate returned data
        R[j,j] = α
        R[j,j+1:end] = 𝐰

        Q.V[j:end, j] = Q₁.v

        Aⱼ = Aⱼ₊₁ # this is the "induction"
        ## END
    end
    Q,R
end

A = randn(600,400)
Q,R = householderqr(A)
@test Q*R ≈ A


# ------ 
# ### Given's Rotations and QR

# An alternative to using reflections to introduce zeros is to use rotations, which
# are called Given's Rotations.
# This is particularly convenient for tridiagonal matrices, where one needs to only
# make one sub-diagonal zero. Here we explore a tridiagonal QR built from rotations
# in a way that the factorisation can be computed in $O(n)$ operations.


# -----


# **Problem 4** This problem explores computing  a QR factorisation of a Tridiagonal matrix in $O(n)$ operations.
# This will introduce entries in the second super-diagonal, hence we will use the `UpperTridiagonal` type
# from Lab 6 (solution copied below). Complete the implementation of `bandedqr`, that only takes $O(n)$ operations,
# using an instance of `Reflections` to represent `Q` and `UpperTriangular` to represent `R`.

import Base: *, size, getindex, setindex!
struct UpperTridiagonal{T} <: AbstractMatrix{T}
    d::Vector{T}   # diagonal entries
    du::Vector{T}  # super-diagonal enries
    du2::Vector{T} # second-super-diagonal entries
end

size(U::UpperTridiagonal) = (length(U.d),length(U.d))

function getindex(U::UpperTridiagonal, k::Int, j::Int)
    d,du,du2 = U.d,U.du,U.du2
    if j - k == 0
        d[j]
    elseif j - k == 1
        du[k]
    elseif j - k == 2
        du2[k]
    else
        0
    end
end

function setindex!(U::UpperTridiagonal, v, k::Int, j::Int)
    d,du,du2 = U.d,U.du,U.du2
    if j > k+2
        error("Cannot modify off-band")
    end
    if j - k == 0
        d[k] = v
    elseif j - k == 1
        du[k] = v
    elseif j - k == 2
        du2[k] = v
    else
        error("Cannot modify off-band")
    end
    U # by convention we return the matrix
end


function bandedqr(A::Tridiagonal)
    n = size(A, 1)
    Q = Rotations(zeros(n - 1)) # Assume Float64
    R = UpperTridiagonal(zeros(n), zeros(n - 1), zeros(n - 2))

    ## TODO: Populate Q and R by looping through the columns of A.

    ## SOLUTION
    ## In what follows we use both A and R simultaneously, where R
    ## represents the upper-triangular part of the modified A, up to rows j.
    ## At each stage we apply the rotation Q_j that introduces a zero
    ## in the (j+1,j) entry of A. 

    ## To begin with we haven't applied Q so R contains just the first row of A.

    R[1, 1:2] = A[1, 1:2]
        
    for j = 1:n-1
        ## We determine the angle of rotation. Note R contains the updated entries
        ## on or above the diagonal whilst A containss the unmodified entries below the diagonal.
        x_1 = [R[j, j], A[j+1, j]]
        Q.θ[j] = atan(x_1[2], x_1[1])
        Q_1 = Rotation(-Q.θ[j]) # rotate in opposite direction 
        ## Q satisfies Q*[ R[j,j],A[j+1,j]] = [sqrt(R[j,j]^2+A[j+1,j]^2), 0]
        ## Or we can just write this as:
        R[j, j] = (Q_1 * x_1)[1]

        ## We now apply this to the rest of the columns to update R.
        x_2 = [R[j, j+1]; A[j+1, j+1]]
        R[j:j+1, j+1] = Q_1 * x_2

        if j < n - 1 # need to avoid going beyond the dimension
            ## We use the fact that A is banded so there is a 0.
            x_3 = [0, A[j+1, j+2]]
            R[j:j+1, j+2] = Q_1 * x_3
        end
    end
    ## END
    Q, R
end

A = Tridiagonal([1, 2, 3, 4], [1, 2, 3, 4, 5], [1, 2, 3, 4])
Q, R = bandedqr(A)
@test Q*R ≈ A



# ### III.6.3 QR and least squares

# When we type `A \ b` with a rectangular matrix `A` it is
# solving a least squares system, and behind the scenes it is using a QR factorisation.
# We can see this via the inbulit `qr` function 

A = randn(200,100)
b = randn(200)

Q,R̂ = qr(A)

# Here `Q` is a special type representing an orthogonal matrix.
# `R̂` is an `UpperTriangular`, that is, we only store the upper triangular
# entries of `R` (which is the same as the reduced QR factorisation). 
# Thus to solve a least squares problem we need to drop the extra entries as
# follows:

c = Q'b # invert Q
c̃ = c[1:size(R̂,1)] # drop extra entries
@test A \ b ≈ R̂ \ c̃

# **Problem 5** Complete the function `leastsquares(A, b)` that uses your
# `householderqr` function to solve a least squares problem.

function leastsquares(A, b)
    ## TODO: use householderqr to solve a least squares problem.
    ## SOLUTION
    m,n = size(A)
    Q, R = householderqr(A)
    UpperTriangular(R[1:n,1:n])\(Q'b)[1:n]
    ## END
end

@test A\b ≈ leastsquares(A,b)


# ## III.4 Polynomial Interpolation and Regression

# We now explore the practical usage of polynomial interpolation and regression.
# In particular we will see that polynomial interpolation may fail as the number
# of points becomes large. 

# ### III.4.1 Polynomial Interpolation

# A quick-and-dirty way to to do interpolation is to invert the Vandermonde matrix.
# That is, for
# $$
# p(x) = ∑_{k = 0}^{n-1} c_k x^k
# $$
# and $x_1, …, x_n ∈ ℝ$, we choose $c_k$ so that $p(x_j) = f(x_j)$ for
# $j = 1, …, n$. We do so by creating the square Vandermonde matrix 
# $$
# V := \begin{bmatrix} 1 & x_1 & ⋯ & x_1^{n-1} \\
#                     ⋮ & ⋮ & ⋱ & ⋮ \\
#                     1 & x_n & ⋯ & x_n^{n-1}
#                     \end{bmatrix}.
# $$
# If the function samples are
# $$
#  𝐟 = \begin{bmatrix} f(x_1) \\ ⋮ \\ f(x_n) \end{bmatrix}
# $$
# then the coefficients of the interpolatory polynomial
# $$
#       𝐜 = \begin{bmatrix}
#           c_0 \\ ⋮ \\ c_{n-1} \end{bmatrix} 
# $$
# must satisfy $V 𝐜 = 𝐟$.  Thus inverting the Vandermonde matrix tells us the coefficients.

# Here we see an example of this using `n` evenly spaced points:

f = x -> cos(10x)
n = 5
𝐱 = range(0, 1; length=n) # evenly spaced points (BAD for interpolation)
V =  [𝐱[j]^k for j = 1:n, k = 0:n-1] # Vandermonde matrix, also can be written as x .^ (0:n)'
𝐟 = f.(𝐱) # evaluate f at x[k], equivalent to [f(x[k]) for k = 1:n]
𝐜 = V \ 𝐟 # invert the Vandermonde matrix and determine the coefficients
p = x -> dot(𝐜, x .^ (0:n-1)) # take a dot product with monomials x .^ 0:n-1 == [x^j for j=0:n-1]
@test p.(𝐱) ≈ V * 𝐜 # evaluating the polynomial on x is the same as applying V


𝐠 = range(0,1; length=1000) # plotting grid, sample a lot more than interpolation points

## To evaluate a polynomial on the plotting grid its faster to create the rectangular Vandermonde matrix associated with that grid:
V_g = [𝐠[j]^k for j = 1:length(𝐠), k = 0:n-1]

plot(𝐠, f.(𝐠); label="function")
plot!(𝐠, V_g*𝐜; label="interpolation")
scatter!(𝐱, f.(𝐱); label="samples")


# Whether an interpolation is actually close to a function is a subtle question,
# involving properties of the function, distribution of the sample points $x_1,…,x_n$,
# and round-off error.
# A classic example is:
# $$
#   f_M(x) = {1 \over M x^2 + 1}
# $$
# where the choice of $M$ can dictate whether interpolation at evenly spaced points converges.

# -------

# **Problem 6(a)** Interpolate $1/(4x^2+1)$ and $1/(25x^2 + 1)$ at an evenly spaced grid of $n$
# points, plotting the solution at a grid of $1000$ points. For $n = 50$ does your interpolation match
# the true function?  Does increasing $n$ to 400 improve the accuracy? How about using `BigFloat`?
# Hint: make sure to make your `range` be `BigFloat` valued, e.g., `range(big(-1), big(1); length=n)`.

## TODO: interpolate 1/(10x^2 + 1) and 1/(25x^2 + 1) at $n$ evenly spaced points, plotting both solutions evaluated at
## the plotting grid with 1000 points, for $n = 50$ and $400$.

## SOLUTION

n = 50
𝐱 = range(-1, 1; length=n)
𝐠 = range(-1, 1; length=1000) # plotting grid

V = 𝐱 .^ (0:n-1)'
V_g = 𝐠 .^ (0:n-1)'

f_4 = x -> 1/(4x^2 + 1)
𝐜_4 = V \ f_4.(𝐱)
f_25 = x -> 1/(25x^2 + 1)
𝐜_25 = V \ f_25.(𝐱)

plot(𝐠, V_g*𝐜_4; ylims=(-1,1))
plot!(𝐠, V_g*𝐜_25)
## We see large errors near ±1 for both examples. 


n = 400
𝐱 = range(-1, 1; length=n)

V = 𝐱 .^ (0:n-1)'
V_g = 𝐠 .^ (0:n-1)'
f_4 = x -> 1/(4x^2 + 1)
𝐜_4 = V \ f_4.(𝐱)
f_25 = x -> 1/(25x^2 + 1)
𝐜_25 = V \ f_25.(𝐱)

plot(𝐠, V_g*𝐜_4; ylims=(-1,1))
plot!(𝐠, V_g*𝐜_25)
##  M = 4 appears to converge whilst M = 25 breaks down.

## Now do big float
n = 400
𝐱 = range(big(-1), 1; length=n)
𝐠 = range(big(-1), 1; length=1000) # plotting grid

V = 𝐱 .^ (0:n-1)'
V_g = 𝐠 .^ (0:n-1)'

f_4 = x -> 1/(4x^2 + 1)
𝐜_4 = V \ f_4.(𝐱)
f_25 = x -> 1/(25x^2 + 1)
𝐜_25 = V \ f_25.(𝐱)

plot(𝐠, V_g*𝐜_4; ylims=(-1,1))
plot!(𝐠, V_g*𝐜_25)
## With M = 4 it looks like it now is converging. This suggests the issue before was numerical error.
## For M = 25 the solution is even less accurate, which suggests the issue is a lack of mathematical
## convergence.

## END

# ------

# ### III.4.2 Polynomial regression

# To overcome issues with interpolation we will instead use regression: use more points than
# the degree of the polynomial. As an example, suppose we want to fit noisy data by a quadratic
# $$
# p(x) = c₀ + c₁ x + c₂ x^2
# $$
# That is, we want to choose $c₀,c₁,c₂$ at data samples $x_1, …, x_m$ so that the following is true:
# $$
# c₀ + c₁ x_j + c₂ x_j^2 ≈ f_j
# $$
# where $f_j$ are given by data. We can reinterpret this as a least squares problem: minimise the norm
# $$
# \left\| \begin{bmatrix} 1 & x_1 & x_1^2 \\ ⋮ & ⋮ & ⋮ \\ 1 & x_m & x_m^2 \end{bmatrix}
# \begin{bmatrix} p₀ \\ p₁ \\ p₂ \end{bmatrix} - \begin{bmatrix} f_1 \\ ⋮ \\ f_m \end{bmatrix} \right \|
# $$
# When a matrix is rectangular `\` solves a least squares problem for us:

m,n = 100,3

𝐱 = range(0,1; length=m) # 100 points
𝐟 = 2 .+ 𝐱 .+ 2𝐱.^2 .+ 0.1 .* randn.() # Noisy quadratic samples, built with broadcast notation.

V = 𝐱 .^ (0:2)'  # 100 x 3 Vandermonde matrix, equivalent to [ones(m) x x.^2]

𝐜 = V \ 𝐟 # coefficients are, very roughly, [2,1,2]

# We can visualise the fit:

𝐠 =range(0, 1; length=1000)

p = x -> 𝐜[1] + 𝐜[2]x + 𝐜[3]x^2

scatter(𝐱, 𝐟; label="samples", legend=:bottomright)
plot!(𝐠, p.(𝐠); label="quadratic")

# -----

# **Problem 6(b)** Repeat the previous problem but now using _least squares_: instead of interpolating,
# use least squares on a large grid: choose the coefficients of a degree $(n-1)$ polynomial so that
# $$
#     \left\| \begin{bmatrix} p(x_1) \\ ⋮ \\ p(x_m) \end{bmatrix} - \begin{bmatrix} f(x_1) \\ ⋮ \\ f(x_m) \end{bmatrix} \right \|.
# $$
# is minimised, where $n = 50$ and $m = 500$. 
# Does this improve the accuracy near the endpoints? Do you think convergence for a least squares approximation
# is dictated by the radius of convergence of the corresponding Taylor series?
# Hint: use the rectangular Vandermonde matrix to setup the Least squares system.

## TODO: approximate 1/(10x^2 + 1) and 1/(25x^2 + 1) using a least squares system where the 

## SOLUTION
n = 50 # use basis [1,x,…,x^(49)]
𝐱 = range(-1, 1; length=500) # least squares grid
𝐠 = range(-1, 1; length=2000) # plotting grid

V = 𝐱 .^ (0:n-1)'
V_g = 𝐠 .^ (0:n-1)'
f_4 = x -> 1/(4x^2 + 1)
𝐜_4 = V \ f_4.(𝐱)
f_25 = x -> 1/(25x^2 + 1)
𝐜_25 = V \ f_25.(𝐱)

plot(𝐠, V_g*𝐜_4; ylims=(-1,1))
plot!(𝐠, V_g*𝐜_25)

## Yes, now both approximations appear to be converging.
## This is despite the radius of convergence of both functions being
## smaller than the interval of interpolation.

## END
