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