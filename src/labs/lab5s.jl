# # MATH50003 (2024–25)
# # Lab 5: III.2 LU Factorisation and III.3 Cholesky Factorisation

# In this lab we explore using LU, PLU and Cholesky factorisations, and
# implement algorithms for computing a Cholesky factorisation. We explore
# stability properties of these different factorisations, and see that the
# Cholesky factorisation is a robust way of determining if a matrix is symmetric
# postive definite. 



# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 1. LU and PLU factorisations
# 2. Cholesky and reverse Cholesky factorisations, including for banded matrices.

#
# Coding knowledge:
#
# 1. Using the `lu` and `cholesky` functions.


# We load the following packages:

using LinearAlgebra, Plots, Test

# ### III.2 LU and PLU Factorisations

# LU, PLU and Cholesky factorisations are closely related
# matrix factorisations that reduce a square matrix to a product of
# lower and upper triangular matrices, possibly with a permutation matrix.
# We will only focus on the practical usage of LU and PLU, without digging into the
# details of implementation. 

# ### III.2.1 LU Factorisation


# If $A ∈ 𝔽^{n × n}$ is a square matrix where $𝔽$ is a field ($ℝ$ or $ℂ$)
# then we can sometimes find lower and upper triangular matrices $L,U ∈ 𝔽^{n × n}$ such that 
# $$
# A = LU.
# $$
# This is equivalent to Gaussian elimination but we will only focus on practical usage in this lab.
# This factorisation can be computed using the `lu` function, but as the default is a PLU factorisation we add a flag
# telling it not to use pivoting/permutations:

A = [1.0 1 1;
     2   4 8;
     1   4 9]

L,U = lu(A, NoPivot()) # NoPivot is a flag that tells lu to not use permutations

# This matches what we derived by hand in the notes and indeed:

@test A ≈ L*U

# We can use an LU factorisation to reduce solving a linear system to inverting triangular matrices:

b = randn(3)
c = L \ b # computed using forward elimination (even though L is a Matrix, \ detects it is lower triangular)
x = U \ c # computed using back substitution
@test A \ b ≈ x

# If a matrix has a zero on a pivot we know by equivalence to Gaussian elimination that an LU factorisation
# does not exist:

A[1,1] = 0
lu(A, NoPivot()) # throws an error

# But even if it has a very small but non-zero entry we get huge errors:

A[1,1] = 1E-14
L,U = lu(A, NoPivot()) # Succeeds but suddenly U is on order of 2E14!
#
norm(A \ b - U\(L\b)) # Very large error! A \ b uses pivoting now.

# **WARNING** The parantheses are important: algebra is left-associative so had we written `U\L\b` this would have been interpreted as
# `(U\L) \ b` which would have meant `inv(inv(U)*L)*b == L \ (U*b)`.

# -----

# **Problem 1** For `A` defined above, consider setting  `A[1,1] = ε` for `ε = 10.0 ^ (-k)` for `k = 0,…,14`
# with the right-hand side `b = [1,2,3]`.
# Plot, scaling the $y$-axis logarithmically, the growth rate in the error of using LU compared to `\`.
# Make a conjecture on how the error grows as $k → ∞$.
# Hint: you can either allocate a vector of errors that is populated in a for-loop or write a simple comprehension.

## TODO: Do a log-log plot for A with its 1,1 entry set to different ε and guess the growth rate.
## SOLUTION

A = [1.0 1 1;
     2   4 8;
     1   4 9]

b = [1,2,3]


n = 15
errs = zeros(n)
for k = 1:n
    A[1,1] = 10.0 ^ (1-k)
    L,U = lu(A, NoPivot())
    errs[k] = norm(A\b - U \ (L \ b))
end

nanabs(x) = x == 0 ? NaN : abs(x)
scatter(0:n-1, nanabs.(errs); yscale=:log10, xticks = 0:15, yticks= 10.0 .^ (-16:0), legend=:bottomright)
## The error grows exponentially, roughly like 10^k

## END



# -----

# ## III.2.2 PLU Factorisation

# In general it is necessary to use pivoting, a feature you have seen
# in Gaussian elimination but as Problem 1 demonstrates we need to do so even if we do not encounter
# a zero. This corresponds to a factorisation of the form
# $$
#  A = P^⊤LU
# $$
# where $P$ is a permutation matrix, $L$ is lower triangular and $U$ is upper triangular.
# We compute this as follows, printing out the permutation:

A = [0.1 1 1;
     2   4 8;
     1   4 9]

L,U,σ = lu(A)
σ

# The permutation matrix is encoded as a vector $σ$. More precisely, we have
# $$
#     P 𝐯 = 𝐯[σ]
# $$
# Thus we can solve a linear system by  first permuting the entries of the right-hand side:

b = [10,11,12]
b̃ = b[σ] # permute the entries to [b[2],b[3],b[1]]

# And then inverting $L$ and $U$ as before:
c = L \ b̃ # invert L with forward substitution
x = U \ c # invert U with back substitution

@test x == A \ b # \ also use PLU to do the solve so these exactly equal

# Note in the following problems we will see that PLU is _usually_ stable but not always.
# Fortunately the set of matrices where it fails to be accurate has extremely small measure.
# The big _open problem_ in numerical linear algebra is to turn this observation into a precise statement.

    

# -----


# **Problem 2(a)** Complete the function  `badmatrix(n)` that returns the following $ℤ^{n × n}$ matrix:
# $$
#   B_n := \begin{bmatrix}
#       1      &&&& 1  \\
#       -1 & 1       &&& 1   \\
#       ⋮ & ⋱ & ⋱   && ⋮    \\
#       -1 & ⋯ & -1 & 1 & 1 \\
#       -1 & ⋯ & -1 & -1 & 1
#   \end{bmatrix}
# $$
# That is: all entries below the diagonal are $-1$ whilst the diagonal and last column are $1$.
#

function badmatrix(n)
    ## TODO: make the "bad matrix" with `Int` entries defined above and return it
    ## SOLUTION
    A = zeros(Int, n, n)
    for k = 1:n
        A[k,n] = 1
    end
    for j = 1:n-1
        A[j,j] = 1
        for k = j+1:n
            A[k,j] = -1
        end
    end
    A
    ## END
end

@test badmatrix(3) isa Matrix{Int}
@test badmatrix(3) == [1 0 1; -1 1 1; -1 -1 1]


# **Problem 2(b)** Does `lu` use pivoting with `badmatrix(n)`? Does it use
# pivoting with a small random perturbation (created via `randn(n,n)`)?

## TODO: Use `lu` on `badmatrix(n)` and a small perturbation to determine if it
## is using pivoting.

## SOLUTION
lu(badmatrix(5)).p # == 1:5, that is no pivoting has occurred
lu(badmatrix(5) + eps()randn(5,5)).p # ≠ 1:5, we have pivoted
## END


# **Problem 2(c)** We can test the accuracy of a method for inverting a matrix
# by applying the matrix and seeing how different it was from the input, eg.
# computing `norm(A*(A\b) - b)`. This would be zero if everything was done with
# exact arithmetic. Plot the norm of this error for `b = randn(n)` for `bandmatrix(n)` 
# and `badmatrix(n) + 1E-15*randn(n,n)` for `n = 25, 50, 75, 100` and
# compare the observed differences in accuracy of PLU.

## TODO: plot the error norm(A*(A\b) - b) for the perturbed and unperturbed badmatrix(n).
## What do you observe?
## SOLUTION
baderrs = zeros(4)
perterrs = zeros(4)

for k = 1:4
    n = 25k
    b = randn(n)
    A = badmatrix(n)
    Ã = A + 1E-15*randn(n,n)
    baderrs[k] = norm(A * (A \ b) - b)
    perterrs[k] = norm(Ã * (Ã \ b) - b)
end

plot(25:25:100, baderrs; label="bad", yscale=:log10)
plot!(25:25:100, perterrs; label="perturbed")

## The perturbation errors stay very small, whilst the unperturbed
## errors blow up.
## END


# -----
    
# ## III.3 Cholesky factorisation

# The Cholesky factorisation is a special case of LU factorisation for the case
# when a matrix is symmetric positive definite (SPD). Hidden in the proof that a Cholesky factorisation
# exists if and only if the matrix is SPD is a simple algorithm for computing it:


function mycholesky(A)
    T = eltype(A)
    m,n = size(A)
    if n ≠ m
        error("Matrix must be square")
    end
    if A ≠ A'
        error("Matrix must be symmetric")
    end

    L = LowerTriangular(zeros(T,n,n)) # a lower triangular matrix that at the end will satisfy L'L
    Aⱼ = copy(A)
    for j = 1:n
        α,𝐯 = Aⱼ[1,1],Aⱼ[2:end,1]
        if α ≤ 0
            error("Matrix is not SPD") # this error would be a proof that the matrix is not SPD, if done rigorously
        end 
        L[j,j] = sqrt(α)
        L[j+1:end,j] = 𝐯/sqrt(α)

        ## induction part
        K = Aⱼ[2:end,2:end] # drop first row and column of A
        Aⱼ = K - 𝐯*𝐯'/α
    end
    L
end

A = Symmetric(rand(100,100) + 100I) # Symmetric takes in a matrix and produces a symmetric version using the upper-triangular part.
L = mycholesky(A)
@test A ≈ L*L'



# With exact arithmetic algorithm succeeds if and only if $A$ is symmetric positive definite.
# With floating point errors this is not necessarily the case. (One could run it with interval arithmetic
# but that would only prove a matrix is SPD if the algorithm succeeded, failure could be caused by
# rounding.)

# In practice one would normally use the inbuilt `cholesky` function as follows:

L̃ = cholesky(A).L
@test L̃ ≈ L # our implementation matches (approximately) the high-performance implementation.

# In the following problem we consider a Cholesky factorisation for tridiagonal matrices. Since we are assuming the
# matrix is symmetric, we will use a special type `SymTridiagonal` that captures the symmetry.
# In particular, `SymTridiagonal(dv, ev) == Tridiagonal(ev, dv, ev)`.


# -----

# **Problem 3** Use `mycholesky` or `cholesky` to deduce if the following matrices are SPD.
# $$
# \begin{bmatrix} 1 & -1  \\
# -1 & 3
# \end{bmatrix}, \begin{bmatrix} 1 & 2 & 2  \\
# 2 & 1 & 2\\
# 2 & 2 & 1
# \end{bmatrix}, \begin{bmatrix} 3 & 2 & 1  \\
# 2 & 4 & 2\\
# 1 & 2 & 5
# \end{bmatrix}, 
# \begin{bmatrix} 4 & 2 & 2 & 1  \\
# 2 & 4 & 2 & 2\\
# 2 & 2 & 4 & 2 \\
# 1 & 2 & 2 & 4
# \end{bmatrix}
# $$

## TODO: Check if you got PS5 Q1 correct using a computer to do the Cholesky factorisation.
## SOLUTION
cholesky([1 -1; -1 3]) # succeeds so is SPD
cholesky([1 2 2; 2 1 2; 2 2 1]) # throws an error so not SPD
cholesky([3 2 1; 2 4 2; 1 2 5]) # succeeds so is SPD
cholesky([4 2 2 1; 2 4 2 2; 2 2 4 2; 1 2 2 4]) # succeeds so is SPD
## END


# **Problem 4** Complete the following
# implementation of `mycholesky` to return a `Bidiagonal` cholesky factor in $O(n)$ operations.


## return a Bidiagonal L such that L'L == A (up to machine precision)
## You are allowed to change A
function mycholesky(A::SymTridiagonal)
    d = A.dv # diagonal entries of A
    u = A.ev # sub/super-diagonal entries of A
    T = float(eltype(A)) # return type, make float in case A has Ints
    n = length(d)
    ld = zeros(T, n) # diagonal entries of L
    ll = zeros(T, n-1) # sub-diagonal entries of L

    ## TODO: populate the diagonal entries ld and the sub-diagonal entries ll
    ## of L so that L*L' ≈ A
    ## SOLUTION
    ld[1] = sqrt(d[1])
    for k = 1:n-1
        ll[k] = u[k]/ld[k]
        ld[k+1] = sqrt(d[k+1]-ll[k]^2)
    end
    ## END

    Bidiagonal(ld, ll, :L)
end

n = 1000
A = SymTridiagonal(2*ones(n),-ones(n-1))
L = mycholesky(A)
@test L isa Bidiagonal
@test L*L' ≈ A
#----



