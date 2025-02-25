# # MATH50003 (2024–25)
# # Lab 8: V.1 Fourier Expansions and V.2 Discrete Fourier Transform

# 
# **Learning Outcomes**
#
# Mathematical knowledge:

# 1. 

# Coding knowledge:

# 1. The FFTW.jl package and the `fft` function.
# 2. 


# We first load  packages we need including two new ones, FFTW.jl (for the fast Foourier transform)
# and QuadGK.jl (for black-box numerical integration)

using LinearAlgebra, Plots, FFTW, QuadGK, Test


# ## V.1 Fourier Expansions

# We consider the Fourier expansion
# $$
# f(x) = ∑_{k=-∞}^∞ f̂_k {\rm e}^{{\rm i}kθ}
# $$
# and the closely related Fourier–Taylor expansion
# $$
# f(x) = ∑_{k=0}^∞ f̂_k {\rm e}^{{\rm i}kθ}
# $$
# where
# $$
# f̂_k := {1 \over 2π} ∫_0^{2π} {\rm e}^{-{\rm i}kθ} f(x)  {\rm d}x.
# $$
# We will ultimately use the Trapezium rule to compute $f̂_k$ but for now we will use QuadGK.jl
# to get a high-accuracy approximation. QuadGK.jl implements a black-box algorithm for computing
# integrals via `quadgk`, eg.

σ,ε = quadgk(exp, 0, 1) # integrates exp(x) for x ∈ [0,1]

# This returns a 2-tuple, the first argument is an approximation to the integral

@test σ ≈ exp(1)-1

# whilst the second argument is an estimate for the error, which in this case is spot-on:

abs(σ - (exp(1)-1))

# If we give an extra parameter `atol` we can specify the desired accuracy.
# We can thus use this to compute  approximations to the true Fourier coefficients that are accurate
# to $10^{-12}$ via:

fouriercoefficient(f, k) = quadgk(θ -> f(θ)*exp(-im*k*θ), 0, 2π, atol=1E-12)[1]/(2π)
fouriercoefficient(exp, 0)

# We can use this to approximate a finite truncation of the Fourier series
# $$
# fₙ(θ) = ∑_{k=-m}^m f̂ₖ {\rm e}^{{\rm i}kθ}
# $$

m = 20 # n = 2m+1 Fourier coefficients
f̂ = [fouriercoefficient(exp, k) for k = -m:m]
F = θ -> [exp(im*k*θ) for k = -m:m] # Create the Fourier basis
fₙ = θ -> transpose(F(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, exp.(g); label="f")
plot!(g, real(fₙ.(g)); label="Re fₙ(θ)")
plot!(g, imag(fₙ.(g)); label="Im fₙ(θ)")

# Because of symmetry properties the imaginary part is numerically zero.
# But observe also that the approxiamtion overshoots at $0$ and $2π$, something
# known as _Gibb's phenomenon_.


# **Problem 1(a)** Repeat the above experiment for the non-symmetric truncation
# $$
# fₙ(θ) = ∑_{k=-m}^{m-1} f̂ₖ {\rm e}^{{\rm i}kθ}
# $$
# with $n = 2m$ Fourier coefficients, where $m = 20$, 
# for $f(θ) = \exp(θ)$, $f(θ) = \exp(\cos(θ))$, and $\exp(\exp({\rm i} θ))$. Is the imaginary part still
# numerically zero? # Does the second example appear to converge at $0$ and $2π$? 

m = 20 # n = 2m Fourier coefficients
## TODO: plot the truncated Fourier expansion with coefficients ranging from -m:m-1


# **Problem 1(b)**  Plot the finite Fourier-Taylor expansion 
# $$
# fₙ(θ) = ∑_{k=0}^{ n-1} f̂ₖ {\rm e}^{{\rm i}kθ}
# $$
# where $n = 20$,
# for $f(θ) = \exp(θ)$, $f(θ) = \exp(\cos(θ))$, and $\exp(\exp({\rm i} θ))$.
# For which of these functions does the approximation appear to converge?

n = 20
## TODO: plot the truncated Fourier-Taylor expansion with coefficients ranging from 0:n-1


# **Problem 1(c)** A cosine expansion is a special case of Fourier series of the form
# $$
# f(θ) = ∑_{k=0}^∞ f̌_k \cos k θ.
# $$
# Derive an expression for $f̌_k$ and plot the finite cosine expansion
# $$
# fₙ(θ) = ∑_{k=0}^{ n-1} f̌_k \cos k θ
# $$
# for $f(θ) = \exp(θ)$ and $f(θ) = \exp(\cos(θ))$. You may assume $f$ is real
# and the returned coefficients should be real.
# For which of these functions does the approximation appear to converge?

n = 20
function cosinecoefficient(f, k)
    ## TODO: use quadgk to approximate f̌_k
    
end
@test cosinecoefficient(θ -> exp(cos(θ)), 5) isa Float64
@test cosinecoefficient(θ -> exp(cos(θ)), 5) ≈ 0.0005429263119137845

## TODO: plot the truncated cosine expansion with coefficients ranging from 0:n-1



# ### V.1.1 Convergence of Fourier series



# ## Problem 2
# 



# ## V.2 Discrete Fourier Transform





# The Discrete Cosine Transform and Chebyshev Expansions

# Cosine series and the discrete Cosine transform are close relatives of Fourier series and the discrete Fourier transform.  A _Fourier series_ in triginometric polynomials is an expansion

# $$f(\theta) = \sum_{k=0}^\infty c_k \cos k \theta + \sum_{k=1}^\infty s_k \sin k \theta,$$

# where 

# $$
# \begin{align*}
# c_0 &= {1 \over 2 \pi} \int_0^{2 \pi} f(\theta) \cos k \theta d \theta \cr
# c_k &= {1 \over \pi} \int_0^{2 \pi} f(\theta) \cos k \theta d \theta \cr
# s_k &= {1 \over \pi} \int_0^{2 \pi} f(\theta) \sin k \theta d \theta \cr
# \end{align*}
# $$

# A Cosine series is an expansion like

# $$f(\theta) = \sum_{k=0}^\infty c_k \cos k \theta$$

# That is, it is a Fourier series with only cosine terms, no sine terms.  

# ## Cosine series for even functions

# Even functions satisfy $f(\theta) = f(-\theta)$.  Cosine expansions capture evenness, because each Cosine term is even: $\cos k \theta = \cos(-k\theta)$.


# **Exercise 1(a)** Why is it true that, if $f$ is even, we have

# $$s_k = 0?$$

# Verify this statement for $k=1,2,3$ using `quadgk` for the even function $f(\theta) = e^{\cos \theta}$.



# **Exercise 1(b)** When $f$ is even, the integrals can be reduced to have a period.  Explain why its true that

# \begin{align*}
# c_0 &= {1 \over  \pi} \int_0^{\pi} f(\theta) \cos k \theta d \theta \cr
# c_k &= {2\over \pi} \int_0^{\pi} f(\theta) \cos k \theta d \theta \cr
# \end{align*}

# Verify this statement using `quadgk` for the even function $f(\theta) = e^{\cos \theta}$ for $k=0,1,2,3$.


# **Exercise 1(c)** Use `quadgk` to approximate $c_k$ for $f(\theta) = e^{\cos \theta}$.  Verify

# $$f(\theta) \approx \sum_{k=0}^9 c_k \cos k \theta$$