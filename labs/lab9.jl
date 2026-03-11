# # MATH50003 (2025–26)
# # Lab 9: V.3 Convolutions and VI.1 Orthogonal Polynomials

# **Learning Outcomes**
#
# Mathematical knowledge:

# 1. Using the FFT to compute convolutions.
# 2. The relationship between  convolutions and addition of random variables.
# 3. 

# Coding knowledge:

using FFTW, StatsBase, Plots



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
# C[a]
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
    
end

function *(C::Circulant, x::Vector)
    ## TODO: compute the product of the circulant matrix C with the vector x.
    ## using the FFT to achieve O(n log n) complexity.
    
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
# That is,  consider a random variable $X$ which takes the values $0,…,n-1$ with probabilities $f_1,…,f_n$, which we can do via `sample` from StatsBase.jl:

𝐟 = [1,2,2,3,3]; 𝐟 = 𝐟/sum(𝐟) # normalise the weights
n = length(𝐟)
X = sample(0:n-1, Weights(𝐟), 1000) # 1000 random samples
histogram(X; normalized=:probability, label="X") # plot a histogram, normalised to sum to 1
scatter!(0:n-1, 𝐚; ylims=(0,0.5), label="f") # desired weights

# We can consider the addition with another random variable $Y$ which takes the same values but with probabilities $g_1,…,g_n$:

𝐠 = [4,1,2,2,3]; 𝐠 = 𝐠/sum(𝐠)
Y = sample(0:n-1, Weights(𝐠), 1000)
histogram(Y; normalized=:probability, label="Y")
scatter!(0:n-1, 𝐠; ylims=(0,0.5), label="g")

# We now consider the addition of random variables. In this context, we
# periodicise the value, that is, we consider the random variable
# $$
# X+Y \mod n.
# $$
# We won't prove it, but we can observe numerically that the distribution
# is given precisely by the discrete convolution:

histogram(mod.(X+Y,n); normalized=:probability, label="X+Y mod n")
scatter!(0:n-1, periodicconv(𝐟,𝐠); ylims=(0,0.5), label="f ⋆ g")


# **Problem 2** Consider a discrete sampling of a Gaussian
# $$
# f_k = 
# $$




# ## VI.1 Chebyshev Polynomials

# We now turn to 
# Chebyshev polynomials are orthogonal on $[-1,1]$ with respect to $w(x) = 1/\sqrt{1-x^2}$.
# They actually have an explicit formula as $T_n(x) = \cos n{\rm acos}\, x$. We can plot the first 5:

g = range(-1,1,100) # plot grid
plot(g, cos.((0:4)' .* acos.(g)); label=["T₀" "T₁" "T₂" "T₃" "T₄"])

# These satisfy a simple 3-term recurrence:
# $$
# \begin{align*}
# x T_0(x) &= T_1(x) \\
# x T_n(x) &= T_{n-1}(x)/2 + T_{n+1}(x)/2 
# \end{align*}
# $$
# The recurrence gives us a way of computing $T_n(x)$ without the need for expensive
# trigonometric functions used in the explicit formula:

function chebyshevt(n, x)
    Tₖ₋₁ = 1.0
    if n == 0
        return Tₖ₋₁
    end
    Tₖ = x
    for k = 1:n-1
        Tₖ,Tₖ₋₁ = 2x*Tₖ - Tₖ₋₁, Tₖ
    end
    Tₖ
end

x = 0.1
@test chebyshevt(10, x) ≈ cos(10*acos(x))

# This is particularly useful if you need _all_ of $T_0(x),…,T_n(x)$, eg., to evaluate
# a Chebyshev expansion, as we can save the values in a vector. 


# ------

# **Problem 4**  A Chebyshev expansion is the same as a Cosine expansion of $f(\cos θ)$. Use
# `discretecosinecoefficient` from the previous lab to implement `discretechebyshevcoefficient`
# for computing Chebyshev expansion coefficients.


function discretecosinecoefficient(f, k, n)
    ret = 0.0
    for j = 1:n
        θⱼ = π*(j-1/2)/n
        ret = ret + f(θⱼ)*cos(k*θⱼ)
    end
    if k == 0
        ret/n
    else
        2ret/n
    end
end

function discretechebyshevcoefficient(f, k, n)
    ## TODO: Use discretecosinecoefficient to approximate the Chebyshev coefficients of f
    
end

n = 15
c = [discretechebyshevcoefficient(exp, k, n) for k=0:n-1]
x = 0.1
T = chebyshevt.(0:n-1,x)' # Chebyshev T
@test T*c ≈ exp(x)

# ------

# In the problem sheet we consider the Chebyshev U polynomials defined by
# $$
# U_n(x) =  {\sin(n+1) θ \over \sin θ}
# $$
# which satisfy the recurrence relationship:
# $$
# \begin{align*}
# x U_0(x) &= U_1(x)/2 \\
# x U_n(x) &= {U_{n-1}(x) \over 2} + {U_{n+1}(x) \over 2}.
# \end{align*}
# $$

# **Problem ** Use the recurrence relationship to evaluate $U_n(x)$
# using only floating point arithmetic operaions (i.e., no trigonometry).

function chebyshevu(n, x)
    ## TODO: Implement U_n(x) via the 3-term recurrence relationship
    
end

@test chebyshevu(10, 0.1) ≈ sin(11*acos(0.1))/sin(acos(0.1))

# **Problem 5(b)** In the problem sheet we saw that $T_n'(x) = n U_{n-1}(x)$. Use this to approximate the
# derivative of a function.

function chebyshevdifferentiation(f, n, x)
    ## TODO: Use discretechebyshevcoefficient and the derivative relationship to approximate the derivative of f,
    ## using a Chebyshev expansion of f up to degree n-1
    
end

@test chebyshevdifferentiation(cos, 15, 0.1) ≈ -sin(0.1)

# **Problem 5(c)** For the function $f(x) = {\rm e}^x$ use the previous part
# to approximate the derivative of $f$  at $x = 0$ for increasing choices of `n` up to 50. 
# What do you observe as $n → ∞$ when you plot the errors? 
# Can you achieve higher accuracy than central differences?

## TODO: use chebyshevdifferentiation to approximate the derivative of f and observe the convergence
