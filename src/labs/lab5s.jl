# # MATH50003 (2024–25)
# # Lab 5: III.2 Cholesky Factorisation and III.3 Orthogonal Matrices

# In this lab we explore using LU, PLU and Cholesky factorisations, and
# implement algorithms for computing a Cholesky factorisation. We explore
# stability properties of these different factorisations, and see that the
# Cholesky factorisation is a robust way of determining if a matrix is symmetric
# postive definite. 


# This lab explores orthogonal matrices, including rotations and reflections.
# We will construct special types to capture the structure of these orthogonal operations,
# with the goal of implementing fast matrix*vector and matrix\vector operations.


# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 1. Cholesky and reverse Cholesky factorisations, including for banded matrices.
# 1. Constructing rotation and reflection matrices.
#
# Coding knowledge:
#
# 1. Using the `lu` and `cholesky` functions.


# We load the following packages:

using LinearAlgebra, Plots, Test

# ### III.2 LU and Cholesky Factorisations

# LU, PLU and Cholesky factorisations are closely related
# matrix factorisations that reduce a square matrix to a product of
# lower and upper triangular matrices, possibly with a permutation matrix.
# We will only focus on the practical usage of LU and PLU, without digging into the
# details of implementation. For the Cholesky factorisation we will look at implementation.

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
#     P^⊤ 𝐯 = 𝐯[σ]
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


# **Problem 3(a)** Complete the function  `badmatrix(n)` that returns the following $ℤ^{n × n}$ matrix:
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


# **Problem 3(b)** Does `lu` use pivoting with `badmatrix(n)`? Does it use
# pivoting with a small random perturbation (created via `randn(n,n)`)?

## TODO: Use `lu` on `badmatrix(n)` and a small perturbation to determine if it
## is using pivoting.

## SOLUTION
lu(badmatrix(5)).p # == 1:5, that is no pivoting has occurred
lu(badmatrix(5) + eps()randn(5,5)).p # ≠ 1:5, we have pivoted
## END


# **Problem 3(c)** We can test the accuracy of a method for inverting a matrix
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
    
# ## III.2.3 Cholesky factorisation

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

# **Problem 4** Use `mycholesky` or `cholesky` to deduce if the following matrices are SPD.
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

## TODO: Check if you got PS6 Q1 correct using a computer to do the Cholesky factorisation.
## SOLUTION
cholesky([1 -1; -1 3]) # succeeds so is SPD
cholesky([1 2 2; 2 1 2; 2 2 1]) # throws an error so not SPD
cholesky([3 2 1; 2 4 2; 1 2 5]) # succeeds so is SPD
cholesky([4 2 2 1; 2 4 2 2; 2 2 4 2; 1 2 2 4]) # succeeds so is SPD
## END


# **Problem 5** Complete the following
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



# ## III.5 Orthogonal and Unitary Matrices

# Here we explore representing rotations and reflections, which are
# special types of orthogonal/unitary matrices. 

# ### III.5.1 Rotations

# A (simple) rotation matrix is an element of the special orthogonal group $SO(2)$ and has a matrix representation
# $$
#  \begin{bmatrix} c & -s \\ s & c \end{bmatrix}
# $$
# such that $c^2 + s^2 = 1$. 
# More generally, we can generalise simple rotations on higher dimensional vectors by acting on two indices at a time.
# There are multiple ways of storing a rotation matrix, here we explore the most intuitive (but not the fastest!) way of storing just an angle $θ$
# so that $c = \cos θ$ and $s = \sin θ$.

# We will use a syntax in a struct that forces a field to be a special type. In what follows we define
# the `getindex` by first implementing multiplication, a pattern that will be reused in the problems.



struct Rotation <: AbstractMatrix{Float64}
    θ::Float64 # The ::Float64 means θ can only be a Float64
end

import Base: *, size, getindex

size(Q::Rotation) = (2, 2)

function *(Q::Rotation, x::AbstractVector)
    if length(x) ≠ 2
        error("dimension mismatch")
    end
    θ = Q.θ
    c,s = cos(θ), sin(θ)
    a,b = x # special syntax so that a == x[1] and b == x[2]
    [c*a - s*b, s*a + c*b]
end

function getindex(Q::Rotation, k::Int, j::Int)
    ## We use the overloaded * above as we will follow this pattern later.
    e_k = zeros(2)
    e_j = zeros(2)
    e_k[k] = 1  # will error if k ≠ 1 or 2
    e_j[j] = 1  # will error if j ≠ 1 or 2
    e_k'*(Q*e_j)
end

Q = Rotation(0.1)

# We can test the ability to rotate a vector to the $x$-axis. Here we use the inbuilt `atan(y,x)` function
# to compute the angle of a vector:


x = [-1,-2]
θ = atan(x[2], x[1]) # angle of the vector x
Q = Rotation(-θ) # rotate back
Q * x # first entry is norm(x), second entry is 0


# -----

# **Problem 1** Complete the implementation of `Rotations`, which represents an orthogonal matrix `Q` that is a product
# of rotations of angle `θ[k]`, each acting on the entries `k:k+1`. That is, it returns $Q = Q_1⋯Q_k$ such that
# $$
# Q_k[k:k+1,k:k+1] = 
# \begin{bmatrix}
# \cos θ[k] & -\sin θ[k]\\
# \sin θ[k] & \cos θ[k]
# \end{bmatrix}
# $$
# with all other entries being equivalent to the identity.

struct Rotations <: AbstractMatrix{Float64}
    θ::Vector{Float64} # a vector of angles
end



## we use the number of rotations to deduce the dimensions of the matrix
size(Q::Rotations) = (length(Q.θ)+1, length(Q.θ)+1)

function *(Q::Rotations, x::AbstractVector)
    ## TODO: Apply Q in O(n) operations. You may assume x has Float64 entries.
    ## Hint: you may wish to use copy(x) and only change the relevant entries. 
    ## SOLUTION
    y = copy(x) # copies x to a new Vector 
    θ = Q.θ
    ## Does Q1....Qn x
    for k = length(θ):-1:1
        #below has 4 ops to make the matrix and 12 to do the matrix-vector multiplication,
        #total operations will be 48n = O(n)
        c, s = cos(θ[k]), sin(θ[k])
        y[k:(k+1)] = [c -s; s c] * y[k:(k+1)]
    end
    ## END

    y
end

function getindex(Q::Rotations, k::Int, j::Int)
    ## TODO: Return Q[k,j] in O(n) operations using *.

    ## SOLUTION
    ## recall that A_kj = e_k'*A*e_j for any matrix A
    ## so if we use * above, this will take O(n) operations
    n = size(Q)[1]
    ej = zeros(eltype(Q), n)
    ej[j] = 1
    ## note, must be careful to ensure that ej is a VECTOR
    ## not a MATRIX, otherwise * above will not be used
    Qj = Q * ej
    Qj[k]
    ## END
end

θ = randn(5)
Q = Rotations(θ)
@test Q'Q ≈ I
@test Rotations([π/2, -π/2]) ≈ [0 0 -1; 1 0 0; 0 -1 0]


# ------

# ### III.5.2 Reflections


# We can also construct reflections, defined by a normalised vector $𝐯$ as
# $$
# Q_{𝐯} := I - 2𝐯𝐯^⋆
# $$
# The obvious way is to create a dense vector, eg.

x = randn(5) # vector we want to reflect
v = x/norm(x) # normalise x
Q = I - 2v*v' # a reflection matrix

# Note `I` is a special convenience type that represents the identity matrix for any dimension.

# A special type of reflection is a Householder reflection, which maps a vector to the $x$-axis.
# Using dense matrices we can construct it as follows:


function dense_householderreflection(x)
    y = copy(x)
    if x[1] == 0
        y[1] += norm(x) 
    else # note sign(z) = exp(im*angle(z)) where `angle` is the argument of a complex number
        y[1] += sign(x[1])*norm(x) 
    end
    w = y/norm(y)
    I - 2*w*w'
end


x = randn(3) + im*randn(3)
Q = dense_householderreflection(x)
Q * x # all the entries apart from the first are numerically zero

# A matrix-vector product is $O(n^2)$ operations but we know we can reduce it to $O(n)$.
# Thus we will create a special type to represent the reflection and obtain the better complexity
# multiplication. Because we want the matrix to be real when the entries are real we will use
# a special feature called "templating". Here by adding the `{T}` after the type we allow this to
# be either a `Float64` or `ComplexF64` (or indeed a `BigFloat`). We also do some checking
# to make sure that our defining vector is already normalised. 

struct Reflection{T} <: AbstractMatrix{T}
    v::Vector{T} # T can be either a Float64 or ComplexF64
end

function Reflection(v::Vector)
    T = eltype(v) # find the type of the entries of v
    if !(norm(v) ≈ 1)
        error("input must be normalised")
    end
    Reflection{T}(v) # create an instance of Reflection, specifying the entry type
end


## Implementations of Reflection * AbstractMatrix
## You may wish to use this below to solve Problem 3.
function *(Q::Reflection, X::AbstractMatrix)
    T = promote_type(eltype(Q), eltype(X))
    m,n = size(X)
    ret = zeros(T, m, n)
    for j = 1:n
        ret[:,j] = Q * X[:,j]
    end
    ret
end


# -----

# **Problem 2(a)** Complete the implementation of a type representing an n × n
# reflection that supports `Q[k,j]` in $O(1)$ operations and `*` in $O(n)$ operations.
# The reflection may be complex (that is, $Q ∈ U(n)$ is unitary).

## Represents I - 2v*v'


size(Q::Reflection) = (length(Q.v),length(Q.v))

## getindex(Q, k, j) is synonym for Q[k,j]
function getindex(Q::Reflection, k::Int, j::Int)
    ## TODO: implement Q[k,j] == (I - 2v*v')[k,j] but using O(1) operations.
    ## Hint: the function `conj` gives the complex-conjugate
    ## SOLUTION
    if k == j
        1 - 2Q.v[k]*conj(Q.v[j])
    else
        - 2Q.v[k]*conj(Q.v[j])
    end
    ## END
end
function *(Q::Reflection, x::AbstractVector)
    ## TODO: implement Q*x, equivalent to (I - 2v*v')*x but using only O(n) operations
    ## SOLUTION
    x - 2*Q.v * dot(Q.v,x) # (Q.v'*x) also works instead of dot
    ## END
end

## If your code is correct, these "unit tests" will succeed
n = 10
x = randn(n) + im*randn(n)
v = x/norm(x)
Q = Reflection(v)
@test Q == I-2v*v'
@test Q'Q ≈ I


## We can scale to very large sizes. here we check the reflection property on an 100_000 matrix:
n = 100_000
x = randn(n) + im*randn(n)
v = x/norm(x)
Q = Reflection(v)
@test Q*x ≈ -x


# **Problem 2(b)** Complete the following implementation of a Housholder reflection  so that the
# unit tests pass, using the `Reflection` type created above.
# Here `s == true` means the Householder reflection is sent to the positive axis and `s == false` is the negative axis.
# You may assume `x` has real entries.

function householderreflection(s::Bool, x::AbstractVector)
    ## TODO: return a Reflection corresponding to a Householder reflection
    ## SOLUTION
    y = copy(x) # don't modify x
    if s
        y[1] -= norm(x)
    else
        y[1] += norm(x)
    end
    Reflection(y/norm(y))
    ## END
end

x = randn(5)
Q = householderreflection(true, x)
@test Q isa Reflection
@test Q*x ≈ [norm(x);zeros(eltype(x),length(x)-1)]

Q = householderreflection(false, x)
@test Q isa Reflection
@test Q*x ≈ [-norm(x);zeros(eltype(x),length(x)-1)]



# **Problem 2(c)**
# Complete the definition of `Reflections` which supports a sequence of reflections,
# that is,
# $$
# Q = Q_{𝐯_1} ⋯ Q_{𝐯_m}
# $$
# where the vectors are stored as a matrix $V ∈ ℂ^{n × m}$ whose $j$-th column is $𝐯_j∈ ℂ^n$, and
# $$
# Q_{𝐯_j} = I - 2 𝐯_j 𝐯_j^⋆
# $$
# is a reflection.


struct Reflections{T} <: AbstractMatrix{T}
    V::Matrix{T} # Columns of V are the householder vectors
end

size(Q::Reflections) = (size(Q.V,1), size(Q.V,1))


function *(Q::Reflections, x::AbstractVector)
    ## TODO: Apply Q in O(mn) operations by applying
    ## the reflection corresponding to each column of Q.V to x
    
    ## SOLUTION
    m,n = size(Q.V)
    for j = n:-1:1
        x = Reflection(Q.V[:, j]) * x
    end
    ## END

    x
end


## Implementations of Reflections * AbstractMatrix
## You may wish to use this below to solve Problem 3.
function *(Q::Reflections, X::AbstractMatrix)
    T = promote_type(eltype(Q), eltype(X))
    m,n = size(X)
    ret = zeros(T, m, n)
    for j = 1:n
        ret[:,j] = Q * X[:,j]
    end
    ret
end


function getindex(Q::Reflections, k::Int, j::Int)
    ## TODO: Return Q[k,j] in O(mn) operations (hint: use *)

    ## SOLUTION
    T = eltype(Q.V)
    m,n = size(Q)
    eⱼ = zeros(T, m)
    eⱼ[j] = one(T)
    return (Q*eⱼ)[k]
    ## END
end

import LinearAlgebra: adjoint
function adjoint(Q::Reflections) # called when calling Q'
    ## TODO: return the adjoint as a Reflections
    ## SOLUTION
    Reflections(Q.V[:,end:-1:1])
    ## END
end

Y = randn(5,3)
V = Y * Diagonal([1/norm(Y[:,j]) for j=1:3])
Q = Reflections(V)
@test Q ≈ (I - 2V[:,1]*V[:,1]')*(I - 2V[:,2]*V[:,2]')*(I - 2V[:,3]*V[:,3]')
@test Q' isa Reflections
@test Q' ≈ (I - 2V[:,3]*V[:,3]')*(I - 2V[:,2]*V[:,2]')*(I - 2V[:,1]*V[:,1]')
@test Q'Q ≈ I