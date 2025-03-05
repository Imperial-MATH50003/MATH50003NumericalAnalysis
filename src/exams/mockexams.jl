# # MATH50003 Numerical Analysis (2024–2025) Mock Computer-based Exam

# Instructions:

# 1. You have 15 mins to read the exam beginning when the invigilators instruct. **DO NOT** write or type anything during this time.
# 2. You have 1 hour to complete the exam beginning when the invigilators instruct. You **MUST STOP** typing when the time is complete.
# 3. When finished, save your work and close Visual Studio Code.
# 4. Re-access WISEflow in SchoolYear and click on ‘Upload Paper’ à From Computer.
# 6. You can now complete submission by clicking the green button "Click here to submit".
# 6. If uploading fails please contact an invigilator.

# Instructions for the exam:

# 1. For each problem, replace the `# TODO` to complete the question.
# The unit tests are provided to help you test your answers, but do not guarantee that the answer is correct.
# 2. Problems are marked A/B/C to indicate difficulty ("A" being most difficult).
# 3. All questions are worth 10 marks. Partial credit will be awarded for reasonable attempts or comments outlining a solution even if the tests
# are not passed.
# 3. If you have technical queries please contact an invigilator.
# 4. You may use existing code from the module Github page
# 5. You **MUST NOT** ask for help online or
# communicate with others within or outside the module.
# Failure to follow these rules will be considered academic misconduct.
# 6. **NO USAGE of AI tools** such as ChatGPT or GitHub Co-Pilot.



# You should use the following packages:

using LinearAlgebra, SetRounding, Test

# **WARNING** It may be necessary to restart the kernel if issues arise. Remember to reload the packages
# when you do so.

# -----

# **Problem 1 (C)** Use second-order divided differences
# with an appropriately chosen $h$ to approximate the second derivative of
# $$
# f(x) = \cos(x^2)
# $$
# at $x = 0.1$ with an error less than $10^{-5}$. Note you are not required to choose a "quasi-optimal"
# value for $h$, as long as your choice achieves the specified accuracy.

function fd2(x)
    ## TODO: implement a second-order finite-difference rule 
    ## to approximate f''(x)
    ## for f(x) = cos(x^2)
    ## with step-size h chosen to get sufficient accuracy
    ## SOLUTION
    h = cbrt(eps())
    f = x -> cos(x^2)
    (f(x + h) - 2f(x) + f(x - h)) / h^2 # 5 marks for similar expression
    ## END
end


@test abs(fd2(0.1) + 2*sin(0.1^2) + 4*0.1^2*cos(0.1^2)) ≤ 1E-5

# **Problem 2 (B)** Implement powers of dual numbers to a float $(a+bϵ)^c$ and
# to a dual number $(a+bϵ)^{c+dϵ}$, in a way that is consistent with a "dual-extension",
# e.g. if $f(x) = x^{3/2}$ or $f(x) = x^x$ then we want to define the power function so that
# in both cases $f(a + bϵ) = f(a) + bf'(a)ϵ$.
# Hint: for the second part recall $x^y = \exp(y \log x)$ which reduces the problem
# to single-argument functions where the "dual-extension" is easy to define.

## Represents a + b*ε
struct Dual
    a
    b
end

import Base: ^, *, isapprox
*(x::Dual, y::Dual) = Dual(x.a*y.a, x.a*y.b + x.b*y.a)
isapprox(x::Dual, y::Dual) = x.a ≈ y.a && x.b ≈ y.b # used in tests

function ^(x::Dual, c::Real)
    ## TODO: Implement Dual(a,b)^c returning a Dual whose b-component is consistent
    ## with differentiation.
    ## SOLUTION
    ## 2 marks for returning a Dual
    ## 2 marks if the first component is x.a^c
    Dual(x.a^c, x.b * c * x.a^(c-1))
    ## END
end

@test Dual(1.0,2.0)^0.0 == Dual(1.0, 0.0)
@test Dual(1.0,2.0)^0.5 == Dual(1.0, 1.0)
@test Dual(1.0,2.0)^(-0.5) == Dual(1.0, -1.0)

function ^(x::Dual, y::Dual)
    ## TODO: Implement Dual(a,b)^Dual(c,d), returning a Dual in a way that is consistent with
    ## differentiation: i.e. for the function f(x) = x^x, f(Dual(2,1)) should return
    ## Dual(f(2), f′(2)) where f′(x) denotes the derivative of f.
    ## SOLUTION
    ## We have (a+bε)^(c+dε) == exp((c+dε) * log(a+bε)) == exp((c+dε) * (log(a)+b/a * ε))
    ## == exp(c*log(a) +(d*log(a) + c*b/a)*ε)
    ## == a^c + (d*log(a) + c*b/a)*a^c * ε
    ## 2 marks for returning a Dual
    ## 2 marks if the first component is x.a^y.a
    Dual(x.a^y.a,  x.a^y.a*(y.b*log(x.a) + y.a *x.b/x.a))
    ## END
end


@test Dual(2.0, 1.0) ^ Dual(3.0, 1.0) ≈ Dual(8.0,8*(3/2 + log(2)))

# **Problem 3 (A)** Consider a 2nd order version of a dual number:
# $$
# a + b ϵ_1 + c ϵ_2
# $$
# such that
# $$
# \begin{align*}
# ϵ_1^2 &= ϵ_2, \\
# ϵ_2^2 &= ϵ_1 ϵ_2 =  0.
# \end{align*}
# $$
# Complete the following implementation supporting `+` and `*`. You may
# assume `a,b,c` are `Float64`. Hint: you may need to work out on paper
# how to multiply `(a + b*ϵ_1 + c*ϵ_2)*(d + e*ϵ_1 + f*ϵ_2)` using the
# relationship above.

import Base: *, +, ^
struct Dual2
    a
    b
    c
end

function +(s::Dual2, t::Dual2)
    ## TODO: Implement Dual2(...) + Dual2(...), returning a Dual2
    ## SOLUTION
    Dual2(s.a + t.a, s.b + t.b, s.c + t.c) # 2 marks
    ## END
end

function +(s::Dual2, c::Real)
    ## TODO: Implement Dual2(...) + c, returning a Dual2
    ## SOLUTION
    Dual2(s.a + c, s.b, s.c) # 2 marks
    ## END
end

function *(c::Number, s::Dual2)
    ## TODO: Implement c * Dual2(...), returning a Dual2
    ## SOLUTION
    Dual2(c * s.a, c * s.b, c * s.c) # 2 marks
    ## END
end

function *(s::Dual2, t::Dual2)
    ## TODO: Implement Dual2(...) * Dual2(...), returning a Dual2

    ## SOLUTION
    ## we deduce (s.a + s.b*ϵ_1 + s.c*ϵ_2)*(t.a + t.b*ϵ_1 + t.c*ϵ_2) == 
    ## s.a * t.a + (s.a*t.b + s.b*t.a)*ϵ_1 + (s.a*t.c + s.b*t.b + s.c*t.a)*ϵ_2
    Dual2(s.a * t.a, s.a * t.b + s.b * t.a, s.a * t.c + s.b * t.b + s.c * t.a) # 4 marks
    ## END
end

f = x -> x*x*x + 2x + 1
x = 0.1
@test f(Dual2(x,1.,0.)) == Dual2(f(x), 3x^2+2, 6x / 2)

## This has computed the first and second derivatives as
## f(x) + f'(x)*ϵ_1 + f''(x)/2*ϵ_2 == (x^3 + x) + (3x^2+1)*ϵ_1 + 6x/2*ϵ_2

# **Problem 4 (C)**  Implement the function `issub` that determines whether a `Float16` is a sub-normal number.
# DO NOT use the inbuilt routine `issubnormal`. 


function issub(x::Float16)
    ## TODO: return true if x is a sub-normal float. Otherwise return false
    ## SOLUTION
    σ = 15
    abs(x) < 2.0^(1-σ)
    ## Note other solutions may use strings and bitstring.
    ## 5 marks for using bitstring even if answer is wrong, or for recognising that the number
    ## has to be smaller than a specified amount.
    ## END
end

@test issub(Float16(0))
@test issub(nextfloat(Float16(0)))
@test issub(prevfloat(Float16(0)))
@test !issub(Float16(1))
@test !issub(reinterpret(Float16,0b0000010000000000))
@test issub(reinterpret(Float16,0b0000001111111111))


# **Problem 5 (C)** Complete the following function `divideby3(x)` that
# returns a tuple `a,b` such that `a` is the largest `Float64` less 
# than or equal to `x/3` and `b` is the smallest `Float64` greater than or equal to `x/3`,
# using the `setrounding` function. You may assume the input is a `Float64`.

function divideby3(x)
    ## TODO: return a tuple a,b so that a ≤ x ≤ b where b is either equal to a or the next float
    ## SOLUTION
    ## 4 marks for using setrounding
    a = setrounding(Float64, RoundDown) do
        x/3
    end
    b = setrounding(Float64, RoundUp) do
        x/3
    end
    a,b
    ## END
end

x = 0.1 # arbitary x
a,b = divideby3(x)
@test a ≤ big(x)/3 ≤ b
@test b == nextfloat(a)



# **Problem 6 (B)** Complete the implementation of `LowerTridiagonal` which represents an $n × n$ banded matrix with
# bandwidths $(l,u) = (2,0)$ by storing only its diagonal, sub-diagonal, and second-sub-diagonal as vectors.
# Overload `getindex`, as well as `*` to multiply a lower tridiagonal matrix times a vector in $O(n)$ operations. 

import Base: getindex,  size, *

struct LowerTridiagonal <: AbstractMatrix{Float64}
    d::Vector{Float64}   # diagonal entries of length n
    dl::Vector{Float64}  # sub-diagonal entries of length n-1
    dl2::Vector{Float64} # second-sub-diagonal entries of length n-2
end

size(L::LowerTridiagonal) = (length(L.d),length(L.d))

function getindex(L::LowerTridiagonal, k::Int, j::Int)
    d, dl, dl2 = L.d, L.dl, L.dl2
    ## TODO: return L[k,j].
    ## SOLUTION
    ## 3 marks, give 1 mark if different cases considered
    if k == j
        d[k]
    elseif k == j+1
        dl[j]
    elseif k == j+2
        dl2[j]
    else
        0.0
    end
    ## END
end


function *(L::LowerTridiagonal, x::AbstractVector)
    ## TODO: Return L*x but computed in O(n) operations
    ## SOLUTION
    n,m = size(L)
    b = zeros(n) # returned vector

    ## 7 marks, only give 4 marks if not O(n)
    for j = 1:n, k = j:min(j+2,n) ## 6 marks if k or j vary similar to this.
        b[k] += L[k,j]*x[j]
    end
    b
    ## END
end

n = 10
d, dl, dl2 = randn(n), randn(n-1), randn(n-2)
L = LowerTridiagonal(d, dl, dl2)
@test L == diagm(0 => d, -1 => dl, -2 => dl2)
x = randn(n)
@test L*x ≈ diagm(0 => d, -1 => dl, -2 => dl2)*x


# **Problem 7 (B)** Complete the function `airy(n)` that returns a length-$n+1$ vector
# $$
# \begin{bmatrix}
# u_0 \\
# ⋮ \\
# u_n
# \end{bmatrix}
# $$
# such that $u_k ≈ u(x_k)$ approximates the solution to the equation
# $$
# \begin{align*}
# u(0) &= 1 \\
# u(1) &= 0 \\
# u'' - x u &= 0
# \end{align*}
# $$
# at the point $x_k = k/n$ using the second order finite-difference approximation:
# $$
# u''(x_k) ≈ {u_{k-1} - 2u_k + u_{k+1} \over h^2}
# $$
# for $k = 1, …, n-1$. The solution must use only $O(n)$ operations.


function airy(n)
    ## TODO: return a Vector{Float64} approximating the solution to the ODE
    ## SOLUTION
    ## The answer could also be in terms of a Tridiagonal solve.
    ## If the solve involves a dense matrix than it is 8 marks.
    x = range(0, 1; length = n+1)
    h = step(x)
    ## 4 marks for using Tridiagonal or SymTridiagonal
    ## 2 marks if using Bidiagonal by mistake.
    A = SymTridiagonal(fill(-2 / h^2, n - 1) .- x[2:end-1], fill(1 / h^2, n - 2))
    u = A \ [-1 / h^2; zeros(n - 2)]
    [1; u; 0] ## 2 marks for correct boundary values.
    ## END
end

u = airy(999)
@test length(u) == 1000
@test u[1] == 1
@test u[end] == 0
## this compares agianst the exact formula
@test abs(u[500] - 0.4757167332829094) ≤ 2E-8



# **Problem 8 (C)** Approximate $\exp x$ by a cubic polynomial by minimising
# the least squares error when sampled at $n$ evenly spaced points in $[0,1]$,
# that is, $x_k = (k-1)/(n-1)$ for $k = 1,…,n$,
# returning the coefficients in the monomial basis.


function expfit(n)
    ## TODO: return the coefficients [c_0,c_1,c_2,c_3] of the polynomial
    ## c_0 + c_1*x + c_2*x^2 + c_3*x^3 that minimises the L^2 error at n
    ## evenly spaced samples
    ## SOLUTION
    x = range(0,1; length=n)
    V = x .^ (0:3)' # 4 marks for making a Vandermonde matrix
    V \ exp.(x) # 4 marks for using \, or computing a QR factorisation
    ## END
end

c = expfit(1000)
x = 0.1
@test abs(c[1] + c[2]*x + c[3]*x^2 + c[4]*x^3 - exp(x)) ≤ 1E-3


# **Problem 9.1 (B)** Implement `lowerhouseholderreflection(x)` where `x` is a vector representing $𝐱 ∈ ℝ^n$
#  to return a reflection $Q$  satisfying $Q 𝐱 = -\| 𝐱 \| 𝐞_n$.
#  The function `lowerhouseholderreflection(x)` should return a `Matrix{Float64}`.
# You may assume that `x` is a `Vector{Float64}`.


function lowerhouseholderreflection(x)
    ## TODO: implement the householder reflector defined above
    ## SOLUTION
    ## 10 marks
    y = copy(x)
    y[end] += norm(x) # Only 1 mark if copied unmodified from Lab 7 householderreflection
    w = y/norm(y)
    I - 2w*w'
    ## END
end
x = [1.0,2,3,4]
Q = lowerhouseholderreflection(x)
@test Q*x ≈ [zeros(3); -norm(x)]
@test Q'Q ≈ I
@test Q ≈ Q'

# **Problem 9.2 (A)** Complete the function `ql(A)` that
# returns a QL decomposition, that is, `Q` is an orthogonal
# matrix and `L` is lower triangular satisfying (up to rounding) `A == Q*L`. You may assume that `A`
# is a square `Matrix{Float64}`. Hint: use the previous part to lower triangularise.

function ql(A)
    m,n = size(A)
    m == n || error("not square")
    ## TODO Create Q and L such that Q'Q == I, L is lower triangular, and Q*L ≈ A
    ## SOLUTION
    ## 1 mark for setup
    Aⱼ = A
    Q = Matrix(1.0I, n, n)
    L = zeros(m, n)
    for j = n:-1:1 # 4 marks for recognising the order has to be reversed
        𝐚ⱼ = Aⱼ[:,end] # first columns of Aⱼ, 1 mark for working with this
        Qⱼ = lowerhouseholderreflection(𝐚ⱼ) # 2 marks for calling lowerhouseholderreflection
        QⱼAⱼ = Qⱼ*Aⱼ # multiply Aⱼ by lower Householder reflection
        α,𝐰 = QⱼAⱼ[end,end],QⱼAⱼ[end,1:end-1]

        ## populate returned data, 1 mark
        L[j,j] = α
        L[j,1:j-1] = 𝐰

        ## following is equivalent to Q = Q*[Qⱼ 0; 0 I], 3 marks
        Q[:,1:j] = Q[:,1:j]*Qⱼ

        Aⱼ = QⱼAⱼ[1:end-1,1:end-1] # this is the "induction", 1 mark
    end
    Q,L
    ## END
end

A = [1.0 2 3; 1 4 9; 1 1 1]
Q,L = ql(A)
@test Q'Q ≈ I
@test Q*L ≈ A
@test L ≈ tril(L) # it is acceptable to have small non-zero entries in L

