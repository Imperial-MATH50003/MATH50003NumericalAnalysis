# # MATH50003 (2025–26)
# # Lab 9: V.3 Convolutions and VI.1 Orthogonal Polynomials

# **Learning Outcomes**
#
# Mathematical knowledge:

# 1. Using the FFT to compute convolutions.
# 2. The relationship between  convolutions and addition of random variables.
# 2. Understanding the definition of orthogonal polynomials and their properties.
# 3. Computing orthogonal polynomials using the Gram–Schmidt process.

# Coding knowledge:

using FFTW



# ## V.3 Convolutions

# We saw in lectures that discrete convolutions, circulant matrices, and the FFT are intrinsically linked. In particular, we can recast a discrete convolution defined via
# $$y_k = 𝐚 ⋆ 𝐛 = \sum_{j=0}^{n-1}  a_{k-j} b_j$$
# where view
# $$
# 𝐚 = \begin{bmatrix} a_0 \\ ⋮\\ a_{n-1}\end{bmatrix},\qquad 𝐛 = \begin{bmatrix} b_0 \\ ⋮\\ b_{n-1}\end{bmatrix}
# $$
# and we use the convention that $a_{-j} = a_{n-j}$. Here is a simple implementation using `mod` to handle the wrap-around:

function periodicconv(𝐚, 𝐛)
    n = length(𝐚)
    y = zeros(eltype(𝐚), n)
    for k in 1:n
        for j in 1:n
            y[k] += 𝐚[mod(k-j, n) + 1] * 𝐛[j]
        end
    end
    return y
end

𝐚 = [1,2,3]
𝐛 = [4,5,6]

@test periodicconv(𝐚, 𝐛) ≈ [1*4 + 3*5 + 2*6, 2*4 + 1*5 + 3*6, 3*4 +2*5 + 1*6]


# In the notes we saw that the convolution can be related to the DFT. In particular, we have
# $$
# 𝐚 ⋆ 𝐛 = Q_n^⋆diag(Q_n 𝐚)Q_n𝐛
# $$
# Here we verify this formula:

@test periodicconv(𝐚, 𝐛) ≈ ifft(fft(𝐚) .*  fft(𝐛))

# Note the latter involves complex `Float64` arithmetic, hence its possible that integer results may not be exact.

# We saw that convolutions could also be viewed as matrix-vector products with circulant matrices. In particular, we have
# $$
# L[a]
# $$


# ----
# **Problem 1** Complete the implementation of `Circulant`, which represents a circulant matrix.

struct Circulant <: AbstractMatrix{Float64}
    a::Vector{Float64}
end

import Base: size, getindex, *
size(C::Circulant) = (length(C.a), length(C.a))
function getindex(C::Circulant, k::Int, j::Int)
    ## TODO: return the (k,j) entry of the circulant matrix defined by C.a
    ## You may wish to use the `mod` function to handle the wrap-around.
    ## SOLUTION
    n = length(C.a)
    return C.a[mod(k-j, n) + 1]
    ## END
end

function *(C::Circulant, x::Vector)
    ## TODO: compute the product of the circulant matrix C with the vector x.
    ## using the FFT to achieve O(n log n) complexity.
    ## SOLUTION
    return real(ifft(fft(C.a) .* fft(x))) # we use real since the result should be real, but the FFT may introduce small imaginary parts due to numerical errors.
    ## END
end


𝐚 = [1,2,3]
𝐛 = [4,5,6]
@test Circulant(𝐚) ≈ [1 3 2; 2 1 3; 3 2 1]
@test Circulant(𝐚) * 𝐛 isa Vector{Float64}
@test Circulant(𝐚) * 𝐛 ≈ periodicconv(𝐚, 𝐛) 


# ----

# The convolution over the real line has a natural interpretation in terms of the sum of independent random variables. In particular, if $X$ and $Y$ are independent random variables with probability density functions (PDFs) $f$ and $g$, then the PDF of the sum $Z = X + Y$ is given by the convolution
# $$
# f_Z(x) = (f ⋆ g)(x) = \int_{-\infty}^\infty f(x-y)g(y)dy.
# $$
# This property is also true for discrete random variables, provided that the addition is taken in a periodic sense.






# **Problem 2** 


# ## VI.1 Orthogonal Polynomials

# ## VI.1 General Orthogonal Polynomials

# Orthogonal polynomials are graded polynomials which have the form
# $$
# p_n(x) = k_n x^n + k_n^{(1)} x^{n-1} + ⋯ + k_n^{(n-1)} x + k_n^{(n)}
# $$
# We can store the (currently unknown) coefficients of the orthogonal polynomials as an upper-triangular matrix:
# $$
# R_n = \begin{bmatrix} k_0 & k_1^{(1)} & k_2^{(2)} & ⋯ & k_n^{(n)} \\
#               & k_1 & k_2^{(1)} & ⋯ & k_n^{(n-1)} \\
#                &  & ⋱ & ⋱ & ⋮ \\
#                 & & & k_{n-1} & k_n^{(1)} \\
#                   &&&& k_n
# \end{bmatrix}
# $$
# This can be written as
# $$
# [p_0| …| p_n] = [1| x| …| x^n] R_n
# $$


# If we move $R_n$ to the otherside this can be viewed as a generalisation
# of the QR factorisation:
# $$
# [1| x| …| x^n] = [p_0| …| p_n]  R_n^{-1}
# $$
# And just as Gram–Schmidt can be used to compute a QR factorisation of a matrix,
# we can build monic orthogonal polynomials using a continuous version of Gram–Schmidt:
# $$
#  π_n(x) = x^n - ∑_{k=0}^{n-1} {⟨x^n,π_k ⟩ \over \|π_k\|^2} π_k(x).
# $$
# To deduce $R_n$ from this process we proceed as follows, assuming the inner product is
# $$
# ⟨f,g⟩ := ∫_0^1 f(x) g(x) w(x) {\rm d}x
# $$
# which we approximate with `quadgk`:

function opgramschmidt(w, n)
    R = UpperTriangular(zeros(n,n)) # Connection matrix with monomials
    for j = 1:n
        R[j,j] = 1 # k_j = 1
        for k = 1:j-1
            πₖ = x -> R[1:k,k]'*[x^ℓ for ℓ=0:k-1] # the previously computed OP
            ip = quadgk(x -> x^(j-1) * πₖ(x) * w(x), 0, 1)[1] # ⟨x^n,π_k⟩ 
            nrm = quadgk(x -> πₖ(x)^2 * w(x), 0, 1)[1] # ||π_k||^2. A better version would store this as its repeating the computation for each j
            R[1:k,j] -= ip/nrm * R[1:k,k] # R[1:k,k] gives us the monomial expansion of πₖ
        end
    end
    R
end

# For the special case of $w(x) = 1$ we get:

R = opgramschmidt(x -> 1, 5)

# That is, we have computed the coefficients corresponding to
# $$
# \begin{align*}
# π_0(x) &= 1, \\
# π_1(x) &= x-1/2,\\
# π_2(x) &= x^2 - x + 1/6,\\
# π_3(x) &= x^3 - 3x^2/2 + 3x/5 - 1/20.
# \end{align*}
# $$
# which we computed explicitly in the notes. We can use this to plot the OPs:

g = range(0,1,100) # plot grid
V = g .^ (0:4)' # Vandermonde matrix for monomials
plot(g, V*R; label=["π₀" "π₁" "π₂" "π₃" "π₄"])

# ----

# **Problem 1(a)** Modify `opgramschmidt` to take in the support of the inner product $(α,β)$ and
# not recompute $\|π_k\|^2$ multiple times, and return a tuple containing $R$ and a vector containing $\|π_0\|^2,…,\|π_{n-1}\|^2$.

function opgramschmidt(w, n, α, β)
    R = UpperTriangular(zeros(n,n)) # Connection matrix with monomials
    nrms = zeros(n) # vector of inner products
    ## TODO: Modify the above code to support general weights and not recompute ||π_k||^2
    ## SOLUTION
    for j = 1:n
        R[j,j] = 1 # k_j = 1
        for k = 1:j-1
            πₖ = x -> dot(R[1:k,k],[x^ℓ for ℓ=0:k-1]) # the previously computed OP
            ip = quadgk(x -> x^(j-1) * πₖ(x) * w(x), α, β)[1] # ⟨x^n,π_k⟩ 
            R[1:k,j] -= ip/nrms[k] * R[1:k,k] # R[1:k,k] gives us the monomial expansion of πₖ
        end
        πⱼ = x -> dot(R[1:j,j],[x^ℓ for ℓ=0:j-1]) # the previously computed OP
        nrms[j] =  quadgk(x -> πⱼ(x)^2 * w(x), α, β)[1]
    end
    ## END
    R,nrms
end

R,nrms = opgramschmidt(x -> 1, 3, 0, 1)
@test R ≈ [1 -1/2 1/6;
           0  1   -1;
           0  0    1]
@test nrms ≈ [1,1/12,1/180]

# **Problem 1(b)** Use the new `opgramschmidt` to compute the monic OPs for $\sqrt{1-x^2}$ and $1-x$ on $[-1,1]$
# Do these match the computation from the problem sheet?

## TODO: employ the new opgramschmidt to the two examples from the problem sheet. 
## SOLUTION
opgramschmidt(x -> sqrt(1-x^2), 5, -1, 1)[1]
## Yes it matches 1, x, x^2 - 1/4, x^3 - x/2
opgramschmidt(x -> 1-x, 5, -1, 1)[1]
## Yes it matches 1, x+1/3, x^2 + 2x/5 - 1/5, x^3 + 3x^2/7 - 3x/7 - 3/35
## END

# **Problem 1(c)** By calling `opgramschmidt` implement `orthonormalgramschmidt`
# to return an upper triangular matrix containing the coefficents for
# the orthonormal polynomials expanded in monomials. For the two examples in the previous problem,
# does this match what you derived in the problem sheet?

function orthonormalgramschmidt(w, n, α, β)
    ## TODO: Use opgramschmidt to construct the R for orthonormal polynomials
    ## SOLUTION
    R,nrms = opgramschmidt(w, n, α, β)
    R * Diagonal(nrms .^ (-1/2))
    ## END
end

R = orthonormalgramschmidt(x -> 1, 4, 0, 1)
@test R ≈ [1 -sqrt(3) sqrt(5)   -sqrt(7);
           0 2sqrt(3) -6sqrt(5) 12sqrt(7);
           0 0        6sqrt(5)  -30sqrt(7);
           0 0        0         20sqrt(7)]


# -----

# ### VI.1.1 Three-term recurrence

# Orthogonal polynomials all have three-term recurrences:
# $$
# \begin{align*}
# x p_0(x) &= a_0 p_0(x) + b_0 p_1(x) \\
# x p_n(x) &= c_{n-1} p_{n-1}(x) + a_n p_n(x) + b_n p_{n+1}(x).
# \end{align*}
# $$
# We can use this fact to simplify the Gram-Schmidt process: we do not need
# to orthogonalise against all $pₖ$ if instead of beginning with $x^k$ we begin with
# $x pₖ$. Here we modify the code in `stieltjes` as this is called the _Stieltjes procedure_
# (or sometimes _Lanczos_ which is the name of the corresponding finite-dimensional algorithm): 

function stieltjes(w, n)
    if n < 2
        error("Only works for large enough n")
    end
    R = UpperTriangular(zeros(n,n)) # Connection matrix with monomials
    a = zeros(n-1) # aₖ
    c = zeros(n-2) # cₖ
    R[1,1] = 1

    ## xπ_0 = a_0π_0 + π_1
    ## a_0 = <xπ_0,1>/||π_0||^2
    a[1] = quadgk(x -> x*w(x), 0, 1)[1]/quadgk(w, 0, 1)[1]
    R[2,2] = 1 # coefficients of x*π_0
    R[1,2] = -a[1] # coefficients of -a_0*π_0
    for j = 2:n-1
        R[2:j+1,j+1] = R[1:j,j] # coefficients of x*πⱼ
        πⱼ = x -> R[1:j,j]'*[x^ℓ for ℓ=0:j-1] # π_ⱼ(x)
        πⱼ₋₁ = x -> R[1:j-1,j-1]'*[x^ℓ for ℓ=0:j-2] # π_ⱼ₋₁(x)

        ## xπ_j = c_jπ_{j-1} + a_j*π_j + π_{j+1}
        c[j-1] = quadgk(x -> x*πⱼ(x)*πⱼ₋₁(x)*w(x), 0, 1)[1]/quadgk(x -> w(x)πⱼ₋₁(x)^2, 0, 1)[1]
        R[1:j-1,j+1] -= c[j-1]*R[1:j-1,j-1] # coefficients of -c_{j-1}*π_{j-1}
        a[j] = quadgk(x -> x*πⱼ(x)^2*w(x), 0, 1)[1]/quadgk(x -> w(x)πⱼ(x)^2, 0, 1)[1]
        R[1:j,j+1] -= a[j]*R[1:j,j] # coefficients of -a_j*π_j
    end
    R, a, c
end

R, a, c = stieltjes(x -> 1, 5)
@test R ≈ opgramschmidt(x -> 1, 5)

# Knowing the three-term recurrence actually gives us an explicit formula for the the OPs themselves.
# here is a simple example of evaluating monic OPs:

function monicforward(n, a, c, x)
    πₖ₋₁ = 1.0
    if n == 0
        return πₖ₋₁
    end
    πₖ = x - a[1]

    for k = 2:n
        πₖ,πₖ₋₁ = (x-a[k])*πₖ - c[k-1]*πₖ₋₁,πₖ # a little tuple trick to reuse variables! The RHS is evaluated first and then the variables are updated
    end
    πₖ
end

x = 0.1
@test monicforward(3, a, c, x) ≈ x^3 - 3x^2/2 + 3x/5 - 1/20 # matches explicit formula

# ----

# **Problem 2(a)** Implement `orthonormalstieltjes` for computing the orthonormal polynomials,
# to return an upper triangular matrix containing the coefficents for
# the orthonormal polynomials expanded in monomials and the 3-term recurrences as vectors,
# noting that $c_j = b_j$. 

function orthonormalstieltjes(w, n, α, β)
    if n < 2
        error("Only works for large enough n")
    end
    R = UpperTriangular(zeros(n,n)) # Connection matrix with monomials
    a = zeros(n-1) # aₖ
    b = zeros(n-1) # bₖ

    R[1,1] = 1/sqrt(quadgk(w, α, β)[1])

    ## TODO: complete the implementation populating R, a, b and c.
    ## SOLUTION
    ## xq_0 = a_0q_0 + b_0q_1
    ## first build p_1 via
    ## xq_0 = a_0q_0 + p_1
    ## then divide by ||p_1||.
    a[1] = quadgk(x -> x*w(x)*R[1,1]^2, α, β)[1]

    R[2,2] = R[1,1] # coefficients of x*q_0
    R[1,2] = -a[1]*R[1,1] # coefficients of -a_0*q_0
    p₁ = x -> R[1,2] + R[2,2]*x
    b[1] = sqrt(quadgk(x -> w(x)*p₁(x)^2, α, β)[1])
    R[:,2] /= b[1] # divide by ||p₁||

    for j = 2:n-1
        R[2:j+1,j+1] = R[1:j,j] # coefficients of x*qⱼ
        qⱼ = x -> R[1:j,j]'*[x^ℓ for ℓ=0:j-1] # qⱼ(x)
        qⱼ₋₁ = x -> R[1:j-1,j-1]'*[x^ℓ for ℓ=0:j-2] # qⱼ₋₁(x)

        ## xq_j = b_{j-1}q_{j-1} + a_j*q_j + p_{j+1}
        R[1:j-1,j+1] -= b[j-1]*R[1:j-1,j-1] # coefficients of -b_{j-1}*q_{j-1}
        a[j] = quadgk(x -> x*qⱼ(x)^2*w(x), α, β)[1]
        R[1:j,j+1] -= a[j]*R[1:j,j] # coefficients of -a_j*p_j
        pⱼ₊₁ = x -> R[1:j+1,j+1]'*[x^ℓ for ℓ=0:j] # pⱼ₊₁(x)
        b[j] = sqrt(quadgk(x -> w(x)*pⱼ₊₁(x)^2, α, β)[1]) # ||pⱼ₊₁||
        R[:,j+1] /= b[j]  # divide by ||pⱼ₊₁||
    end
    ## END
    R, a, b
end

R, a, b = orthonormalstieltjes(x -> 1, 5, 0, 1)
@test R ≈ orthonormalgramschmidt(x -> 1, 5, 0, 1)
@test a ≈ fill(0.5,4)
@test b[1:3] ≈ [1/(2sqrt(3)), 1/sqrt(15), 3/(2*sqrt(35))]


# **Problem 2(b)** Implement `forward` for supporting general 3-term recurrences.

function forward(n, a, b, c, x, p₀)
    if n == 0
        return p₀
    end

    ## TODO: implement forward recurrence to compute pₙ, given the recurrence coefficients a,b,c, x and
    ## initial polynomial p₀ as a constant
    ## SOLUTION
    pₖ₋₁ = p₀
    pₖ = (x-a[1])*p₀/b[1]

    for k = 2:n
        pₖ,pₖ₋₁ = ((x-a[k])*pₖ - c[k-1]*pₖ₋₁)/b[k],pₖ # a little tuple trick to reuse variables! The RHS is evaluated first and then the variables are updated
    end
    pₖ
    ## END
end

@test forward(3, a, b, b, 0.1, R[1,1]) ≈ sqrt(7) * (20x^3 - 30x^2 + 12x - 1)

# ---

