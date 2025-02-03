# # MATH50003 (2024–25)
# # Lab 4: II.3 Interval Arithmetic and III.1 Structured Matrices

# This lab explores the usage of rounding modes for floating point arithmetic and how they
# can be used to compute _rigorous_ bounds on mathematical constants such as `ℯ`.
# The key idea is using _interval arithmetic_ to compute the Taylor series which is
# combined with a bound on the error caused by truncating a Taylor series.
# As a fun example, we compute the first 1000 digits of `ℯ`, backed up by a rigorous
# computation.

# We also explore the construction of vectors and matrices, in particular those with sparsity structure
# such as triangular, diagonal, bidiagonal and tridiagonal
# which we capture using special types.

# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 2. Extending interval arithmetic operations to non-positive intervals.
# 3. Combining interval arithmetic with Taylor series bounds for rigorous computations.
# 1. Matrix multiplication, back- and forward-substitution.
# 2. Banded matrices and their utilisation for better complexity linear algebra.
#
# Coding knowledge:
#
# 3. The `promote` command for converting multiple variables to be the same type.
# 4. Using `&&` for "and" and `||` for "or".
# 1. Construction of a dense `Vector` or `Matrix` either directly or via comprehensions or broadcasting.
# 2. The `vec`, `transpose`, `zeros`, `ones`, `fill`, `range`, `promote_type`, and `eltype` functions.
# 3. Using `\` to solve linear systems.
# 4. The `Diagonal`, `Bidiagonal`, and `Tridiagonal` types for making banded matrices.
# 5. Implementing a custom matrix type by subtyping `AbstractMatrix` and overloading `size`, `getindex`, and `setindex!`.


# We need the following packages:

using SetRounding, LinearAlgebra, Test


# -----
#
# ## II.3 Interval Arithmetic

# In lectures we introduced the idea of using set arithmetic to compute rigorous bounds on algorithms.
# To implement this in practice
# we will now create a Type to represent an interval $[a,b] = \{x : a ≤ x ≤ b\}$, which we will call `Interval`.
# We need fields for the left endpoint (`a`) and a right endpoint (`b`):

struct Interval # represents the set [a,b]
    a # left endpoint
    b # right endpoint
end

Interval(x) = Interval(x,x) # Support Interval(1) to represent [1,1]

# For example, if we say `X = Interval(1, 2)` this corresponds to the mathematical interval
# $[1, 2]$, and the fields are accessed via `X.a` and `X.b`.
# We will overload `*`, `+`, `-`, `/` to use interval arithmetic. That is, whenever we do arithmetic with
# an instance of `Interval` we want it to use correctly rounded interval variants. 
# We also need to support `one` (a function that creates an interval containing a single point `1`)
# and `in` functions (a function to test if a number is within an interval).
# To overload these functions we need to import them as follows:

import Base: *, +, -, ^, /, one, in

# We overload `in` as follows:

in(x, X::Interval) = X.a ≤ x ≤ X.b

# The function `in` is whats called an "infix" operation (just like `+`, `-`, `*`, and `/`). We can call it
# either as `in(x, X)` or put the `in` in the middle and write `x in X`. This can be seen in the following:

X = Interval(2.0,3.3)
## 2.5 in X is equivalent to in(2.5, X)
## !(3.4 in X) is equivalent to !in(3.4, X)
2.5 in X, !(3.4 in X)


# We can overload `one` as follows to create an interval corresponding to $[1,1]$.
# The `one(T)` function will create the "multiplicative identity"
# for a given type. For example `one(Int)` will return `1`, `one(Float64)` returns `1.0`,
# and `one(String)` returns "" (because `"" * "any string" == "any string"`):

one(Int), one(Int64), one(String)

# We can also just call it on an instance of the type:

one(2), one(2.0), one("any string")

# For an interval the multiplicative identity is the interval whose lower and upper limit are both 1.
# To ensure its the right type we call `one(X.a)` and `one(X.b)`

one(X::Interval) = Interval(one(X.a), one(X.b))

# Thus the following returns an interval whose endpoints are both `1.0`:

one(Interval(2.0,3.3))

# We now want to overload the operations `+`, `/` and `*` so that we can compute the Taylor
# series of `exp`. We begin with `+`. 

function +(X::Interval, Y::Interval)
    ##
end


+(x::Number, Y::Interval) = Interval(x) + Y # Number is a supertype that contains Int, Float64, etc.
+(X::Interval, y::Number) = X + Interval(y)


## following example was the non-associative example but now we have bounds
Interval(1.1) + Interval(1.2) + Interval(1.3)

## note we are actually doing computations on ${\rm fl}^{nearest}(1.1)$, etc.,
## that is, we haven't accounted in the errors from making the constants. 


# We now implement division, checking that our assumptions 
# are satified. Note that `&&` means "and" whilst `||` means "or",
# While `!` changes a `true` to `false` and vice-versa.


function /(X::Interval, n::Int)
    a,b = promote(X.a, X.b)
    T = typeof(a)
    if !(n > 0 && 0 < a ≤ b)
        error("Input doesn't satisfy positivity assumptions")
    end
    α = setrounding(T, RoundDown) do
            a / n
    end
    β = setrounding(T, RoundUp) do
            b / n
    end
    Interval(α, β)
end

Interval(1.0,2.0)/3 # rounds bottom down and top up

# Finally we overload `*` to behave like the operation `⊗`:

function *(X::Interval, Y::Interval)
    a,b,c,d = promote(X.a, X.b, Y.a, Y.b)
    T = typeof(a)
    if !(0 < a ≤ b && 0 < c ≤ d)
        error("Input doesn't satisfy positivity assumptions")
    end
    α = setrounding(T, RoundDown) do
            a * c
    end
    β = setrounding(T, RoundUp) do
            b * d
    end
    Interval(α, β)
end

# Let's also support powers:

function ^(X::Interval, k::Int)
    if k ≤ 0
        error("not supported")
    elseif k == 1
        X
    else
        X * X^(k-1)
    end
end

# We can now compute positive polynomials with interval arithmetic:

X = Interval(1.0)
1 + X + X^2/2 + X^3/6 + X^4/24

# ------


# **Problem 3(a)** Complete the following implementations of `-` to correctly round
# the endpoints in interval negation and subtraction.

import Base: -

function -(X::Interval)
    a,b = promote(X.a, X.b)
    ## TODO: return an interval representing {-x : x in X}
    
end

function -(X::Interval, Y::Interval)
    a,b,c,d = promote(X.a, X.b, Y.a, Y.b)
    T = typeof(a)
    ## TODO: return an interval implementing X ⊖ Y
    
end

@test -Interval(0.1,0.2) == Interval(-0.2, -0.1)
@test Interval(0.1,0.2) - Interval(1.1,1.2) ≡ Interval(-1.1, -0.9)

# **Problem 3(b)** Alter the implementation of `/(X::Interval, n::Int)`
# to support the case where `n < 0` and `*` to remove the restrictions on
# positivity of the endpoints. You may assume the intervals are non-empty.

## TODO: overload / and *, again.



@test Interval(1.1, 1.2) * Interval(2.1, 3.1) ≡ Interval(2.31, 3.72)
@test Interval(-1.2, -1.1) * Interval(2.1, 3.1) ≡ Interval(-3.72, -2.31)
@test Interval(1.1, 1.2) * Interval(-3.1, -2.1) ≡ Interval(-3.72, -2.31)
@test Interval(-1.2, -1.1) * Interval(-3.1, -2.1) ≡ Interval(2.31, 3.72)


@test Interval(1.0,2.0)/3 ≡ Interval(0.3333333333333333, 0.6666666666666667)
@test Interval(1.0,2.0)/(-3) ≡ Interval(-0.6666666666666667, -0.3333333333333333)

@test Interval(-1., 2) * Interval(2,3) ≡ Interval(-3.0, 6.0)
@test Interval(-1., 2) * Interval(-3,5) ≡ Interval(-6.0, 10.0)

# -----

# The following function  computes the first `n+1` terms of the Taylor series of $\exp(x)$:
# $$
# \sum_{k=0}^n {x^k \over k!}
# $$
# We avoid using `factorial` to avoid underflow/overflow.

##


# In the notes we derived a bound assuming $0 ≤ x ≤ 1$
# on the error in Taylor series of the form $|δ_{x,n}| ≤ 3/(n+1)!$.
# Here we incorporate that error to get a rigorous bound.

##

# Here we test that the bounds match our expectations:

##

# We can even use the code with `BigFloat` to compute a rigorous bound on the first
# 1000 digits of `ℯ`:

##

# Our tests show that this has computed more than 1000 digits:

##


# ------
# **Problem 4** Extend the implementation of `exp_bound` for the case when `-2 ≤ x ≤ 2`.

## TODO: re-overload exp_bound but without the restrictions on positivity and adjusting the
## the bound appropriately.



@test exp(big(-2)) in exp_bound(Interval(-2.0), 20)

# **Problem 5(a)** Complete the implementation of a function `sin_t(x,n)` computing the
# first `2n+1` terms of the Taylor series:
# $$
# \sin\ x ≈ ∑_{k=0}^n {(-1)^k x^{2k+1} \over (2k+1)!}
# $$

function sin_t(x, n)
    ret = x
    s = x
    ## TODO: Compute the first 2n+1 terms of the Taylor series of sin, without using the factorial function
    
    ret
end

@test sin_t(1.0, 10) ≈ 0.8414709848078965
@test sin_t(big(1.0), 10) in  sin_t(Interval(1.0), 10)

# **Problem 5(b)** Complete the implementation of a function `sin_bound(x,n)` that
# includes an error bound on the computation. You may assume $0 ≤ x ≤ 1$.

function sin_bound(X::Interval, n)
    a,b = promote(X.a, X.b)
    T = typeof(a)
    ## TODO: complete the implementation to include the error in truncating the Taylor series. 
    
end


S = sin_bound(Interval(1.0), 20)
@test sin(big(1)) in S
@test S.b - S.a ≤ 1E-13 # we want our bounds to be sharp

# -----

# ## III.1 Structured Matrices

# Before discussing structured matrices we give an overview of creating arrays  (vectors and matrices)
# in Julia.

# ### III.1.1 Dense matrices


# One can create arrays in multiple ways. For example, the function `zeros(Int, 10)` creates
# a 10-element `Vector` whose entries are all `zero(Int) == 0`. Or `fill(x, 10)` creates a
# 10-element `Vector` whose entries are all equal to `x`. Or you can use a comprehension:
# for example `[k^2 for k = 1:10]` creates a vector whose entries are `[1^2, 2^2, …, 10^2]`.
# This also works for matrices: `zeros(Int, 10, 5)` creates a 10 × 5 matrix of all zeros,
# and `[k^2 + j for k=1:3, j=1:4]` creates the following:

##

# Note sometimes it is best to create a vector/matrix and populate it. For example, the
# previous matrix could also been constructed as follows:

##

# **Remark** Julia uses 1-based indexing where the first index of a vector/matrix
# is 1. This is standard in all mathematical programming languages (Fortran, Maple, Matlab, Mathematica)
# whereas those designed for computer science use 0-based indexing (C, Python, Rust).



# Be careful: a `Matrix` or `Vector` can only ever contain entries of the right
# type. It will attempt to convert an assignment to the right type but will throw
# an error if not successful:

##


# ------
# **Problem 1(a)** Create a 5×6 matrix whose entries are `Int` which is
# one in all entries. Hint: use a for-loop, `ones`, `fill`, or a comprehension.
## TODO: Create a matrix of ones, 4 different ways


# **Problem 1(b)** Create a 1 × 5 `Matrix{Int}` with entries `A[k,j] = j`. Hint: use a for-loop or a comprehension.

## TODO: Create a 1 × 5  matrix whose entries equal the column, 2 different ways


# -------
# #### Transposes and adjoints

# We can also transpose a matrix `A` via `transpose(A)`
# or compute the adjoint (conjugate-transpose) via `A'` (which is
# equivalent to a transpose when the entries are real).
# This is done _lazily_: they return custom types `Transpose` or
# `Adjoint` that just wrap the input array and reinterpret the entries.
# This is equivalent to
# _row-major_ format, where the next address in memory of `transpose(A)` corresponds to
# moving along the row.
# Here is a simple example:

##

# If we change entries of `A'` it actually changes entries of `A` too since
# they are pointing to the same locations in memory, just interpreting the data differently:

##

# Note vector adjoints/transposes behave differently than 1 × n matrices: they are
# more like row-vectors. For example the following computes the dot product of two vectors:

##

# #### Broadcasting

# _Broadcasting_ is a powerful and convenient way to create matrices or vectors,
# where a function is applied to every entry of a vector or matrix.
# By adding `.` to the end of a function we "broadcast" the function over
# a vector:

##

# Broadcasting has some interesting behaviour for matrices.
# If one dimension of a matrix (or vector) is `1`, it automatically
# repeats the matrix (or vector) to match the size of another example.
# In the following we use broadcasting to pointwise-multiply a column and row
# vector to make a matrix:

##

# Since `size([1,2,3],2) == 1` it repeats the same vector to match the size
# `size([4,5]',2) == 2`. Similarly, `[4,5]'` is repeated 3 times. So the
# above is equivalent to:
##

# Note we can also use matrix broadcasting with our own functions:
##


# #### Ranges

# _Ranges_ are another useful example of vectors, but where the entries are defined "lazily" instead of
# actually created in memory.
# We have already seen that we can represent a range of integers via `a:b`. Note we can
# convert it to a `Vector` as follows:
##

# We can also specify a step:
##

# Finally, the `range` function gives more functionality, for example, we can create 4 evenly
# spaced points between `-1` and `1`:
##

# Note that `Vector` is mutable but a range is not:
##

# Both ranges `Vector` are subtypes of `AbstractVector`, whilst `Matrix` is a subtype of `AbstractMatrix`.


# -----

# **Problem 1(c)** Create a vector of length 5 whose entries are `Float64`
# approximations of `exp(-k)`. Hint: use a for-loop, broadcasting `f.(x)` notation, or a comprehension.
## TODO: Create a vector whose entries are exp(-k), 3 different ways



# ------
# #### Storage of matrices and vectors

# A `Vector` stores its entries consecutively in memory.
# To be perhaps overly technical: a `Vector` contains a "pointer" (an integer)
# to the first memory address and a length. A `Matrix` is also stored consecutively in memory
#  going down column-by-
# column (_column-major_). That is,
##
# Is actually stored equivalently to a length `6` vector `[A[1,1],A[2,1],A[3,1],A[1,2],A[2,2],A[3,2]]`:
##

# which in this case would be stored using `8 * 6 = 48` consecutive bytes.
# Behind the scenes, a matrix is also "pointer" to the location of the first entry alongside two integers
# dictating the row and column sizes.


# Matrix-vector multiplication works as expected because `*` is overloaded:
##


# We can implement our own version for any types that support `*` and `+` but there are
# actually two different ways. The most natural mathematical way is to multiply-by-rows:
##


# But we can also change the order of operations to give an alternative approach that is multiply-by-columns:
##


# Both implementations match _exactly_ for integer inputs:
##


# Either implementation will be $O(mn)$ operations. However, the implementation
# `mul_cols` accesses the entries of `A` going down the column,
# which happens to be _significantly faster_ than `mul_rows`, due to accessing
# memory of `A` in order. We can see this by measuring the time it takes using `@btime`:
##

# Here `ms` means milliseconds (`0.001 = 10^(-3)` seconds) and `μs` means microseconds (`0.000001 = 10^(-6)` seconds).
# On my machine we observe that `mul_cols` is roughly 2–3x faster than `mul_rows`, while the optimised `*` is roughly 5x faster than `mul_cols`.
# The reason why isn't too important for us (accessing memory in order is much faster than jumping around), but the key points are:
# 1. Making fast algorithms is delicate and arguably more of an art than a science.
# 2. We can focus on complexity rather than counting operations as the latter does not tell us speed.
# 3. Use in-built implementations whenever available.


# Note that the rules of floating point arithmetic apply here. Matrix multiplication with floats
# will incur round-off error (the precise details of which are subject to the implementation):
##
# And integer arithmetic will be subject to overflow:
##
# Solving a linear system is done using `\`:

A = [1 2 3;
     1 2 4;
     3 7 8]
b = [10; 11; 12]
A \ b

# Despite the answer being integer-valued,
# here we see that it resorted to using floating point arithmetic,
# incurring rounding error.
# But it is "accurate to (roughly) 16-digits".
# As we shall see, the way solving a linear system works is we first write `A` as a
# product of matrices that are easy to invert, e.g., a product of triangular matrices or a product of an orthogonal
# and triangular matrix.


# ### III.1.2 Triangular Matrices

# Triangular matrices are represented by dense square matrices where the entries below the
# diagonal are ignored:
##

# We can see that `L` is storing all the entries of `A` in a field called `data`:
##

# Similarly we can create an upper triangular matrix by ignoring the entries below the diagonal:
##

# If we know a matrix is triangular we can do matrix-vector multiplication in roughly half
# the number of operations by skipping over the entries we know are zero:
##

x = [10, 11, 12]
## matches built-in * which also exploits the structure:
@test mul_cols(L, x) == L*x


# Moreover, we can easily invert matrices.
# Consider a simple 3×3 example, which can be solved with `\`:
##

# Behind the scenes, `\` is doing forward-substitution.
# We can implement our own version as follows:

## ldiv is our own version of \
##



# ------

# **Problem 3(a)** Complete the following function for upper triangular matrix-vector
# multiplication without ever accessing the zero entries of `L` above the diagonal.
# You may assume all numbers are `Float64`.
# Hint: just copy code for `mul_cols` and modify the for-loop ranges a la the `UpperTriangular`
# case.

function mul_cols(U::UpperTriangular, x)
    n = size(U,1)

    b = zeros(n) # the returned vector, begins of all zeros

    ## TODO: populate b so that U*x ≈ b
    

    b
end

U = UpperTriangular(randn(5,5))
x = randn(5)
@test U*x ≈ mul_cols(U, x)


# **Problem 3(b)** Complete the following function for solving linear systems with
# upper triangular systems by implementing back-substitution. You may assume
# all input and output vectors have `Float64` values.

## ldiv(U, b) is our implementation of U\b
function ldiv(U::UpperTriangular, b)
    n = size(U,1)

    if length(b) != n
        error("The system is not compatible")
    end

    x = zeros(n)  # the solution vector
    ## TODO: populate x with the entries according to back substitution.
    
    x
end


U = UpperTriangular(randn(5,5))
b = randn(5)
@test U\b ≈ ldiv(U, b)


# ## III.1.3 Banded matrices


# Diagonal matrices in Julia are stored as a vector containing the diagonal entries:
##
# It is clear that we can perform diagonal-vector multiplications and solve linear systems involving diagonal matrices efficiently
# (in $O(n)$ operations).


# We can create bidiagonal matrices in Julia by specifying the diagonal and off-diagonal:
##
# Multiplication and solving linear systems with Bidiagonal systems is also $O(n)$ operations, using the standard
# multiplications/back-substitution algorithms but being careful in the loops to only access the non-zero entries.


# Julia has a type `Tridiagonal` for representing a tridiagonal matrix from its sub-diagonal, diagonal, and super-diagonal:
##
# Tridiagonal matrices will come up in solving second-order differential equations and orthogonal polynomials.
# We will later see how linear systems involving tridiagonal matrices can be solved in $O(n)$ operations.



# -----

# **Problem 4(a)** Complete the implementation of `UpperTridiagonal` which represents a banded matrix with
# bandwidths $(l,u) = (0,2)$ by overloading `getindex(U::UpperTridiagonal, k::Int, j::Int)` (which implements `U[k,j]`) and `setindex!(U::UpperTriangular, v, k::Int, j::Int)` (which implements `U[k,j] = v`). Return zero (of the same type as the other entries)
# if we are off the bands.

struct UpperTridiagonal <: AbstractMatrix{Float64}
    d   # diagonal entries stored as a Vector: d[k] == U[k,k]
    du  # super-diagonal enries stored as a Vector: du[k] == U[k,k+1]
    du2 # second-super-diagonal entries stored as a Vector: du2[k] == U[k,k+2]
end

## This uses the notation `<: AbstractMatrix{Float64}`: this tells Julia that our type is in fact a matrix
## whose entries are `Float64`.
## In order for it to behave a matrix we have to overload the function `size` for our type to return
## the dimensions (in this case we just use the length of the diagonal):

import Base: size, getindex, setindex!

size(U::UpperTridiagonal) = (length(U.d),length(U.d))

## Julia still doesn't know what the entries of the matrix are. To do this we need to overload `getindex`.
## We also overload `setindex!` to allow changing the non-zero entries.


## getindex(U, k, j) is another way to write U[k,j].
## This function will therefore be called when we call U[k,j]
function getindex(U::UpperTridiagonal, k::Int, j::Int)
    d,du,du2 = U.d,U.du,U.du2
    ## TODO: return U[k,j]
    
end

## setindex!(U, v, k, j) gets called when we write (U[k,j] = v).
function setindex!(U::UpperTridiagonal, v, k::Int, j::Int)
    d,du,du2 = U.d,U.du,U.du2
    if j > k+2 || j < k
        error("Cannot modify off-band")
    end

    ## TODO: modify d,du,du2 so that U[k,j] == v
    
    U # by convention we return the matrix
end

U = UpperTridiagonal([1.0,2,3,4,5], [1.0,2,3,4], [1.0,2,3])
@test U == [1 1 1 0 0;
            0 2 2 2 0;
            0 0 3 3 3;
            0 0 0 4 4;
            0 0 0 0 5]

U[3,4] = 2
@test U == [1 1 1 0 0;
            0 2 2 2 0;
            0 0 3 2 3;
            0 0 0 4 4;
            0 0 0 0 5]




# **Problem 4(b)** Complete the following implementations of `*` and `\` for `UpperTridiagonal` so that
# they take only $O(n)$ operations. Hint: the function `max(a,b)` returns the larger of `a` or `b`
# and `min(a,b)` returns the smaller. They may help to avoid accessing zeros.

import Base: *, \

function *(U::UpperTridiagonal, x::AbstractVector)
    n = size(U,1)
    b = zeros(n) # the returned vector, assume Float64 values
    ## TODO: populate b so that U*x ≈ b (up to rounding)
    
    b
end

function \(U::UpperTridiagonal, b::AbstractVector)
    n = size(U,1)
    if length(b) != n
        error("The system is not compatible")
    end

    x = zeros(n)  # the solution vector, assume Float64 values
    ## TODO: populate x so that U*x ≈ b
    
    x
end

n = 1_000_000 # under-scores are like commas: so this is a million: 1,000,000
U = UpperTridiagonal(ones(n), fill(0.5,n-1), fill(0.1,n-2))
x = ones(n)
b = [fill(1.6,n-2); 1.5; 1] # exact result
## note following should take much less than a second
@test U*x ≈ b
@test U\b ≈ x
