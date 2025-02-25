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
## SOLUTION
f̂ = [fouriercoefficient(exp, k) for k = -m:m-1]
F = θ -> [exp(im*k*θ) for k = -m:m-1] # Create the Fourier basis
fₙ = θ -> transpose(F(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, exp.(g); label="f")
plot!(g, real(fₙ.(g)); label="Re fₙ(θ)")
plot!(g, imag(fₙ.(g)); label="Im fₙ(θ)")
## The imaginary part is not zero

f = θ -> exp(cos(θ))
f̂ = [fouriercoefficient(f, k) for k = -m:m-1]
F = θ -> [exp(im*k*θ) for k = -m:m-1] # Create the Fourier basis
fₙ = θ -> transpose(F(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, f.(g); label="f")
plot!(g, real(fₙ.(g)); label="Re fₙ(θ)")
plot!(g, imag(fₙ.(g)); label="Im fₙ(θ)")
## The imaginary part is now numerically zero and the approximation is indistinguishable.


f = θ -> exp(exp(im*θ))
f̂ = [fouriercoefficient(f, k) for k = -m:m-1]
F = θ -> [exp(im*k*θ) for k = -m:m-1] # Create the Fourier basis
fₙ = θ -> transpose(F(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, real.(f.(g)); label="Re f")
plot!(g, imag.(f.(g)); label="Im f")
plot!(g, real(fₙ.(g)); label="Re fₙ") 
plot!(g, imag(fₙ.(g)); label="Im fₙ")
## The real and imaginary parts match to high accuracy
## END

# **Problem 1(b)**  Plot the finite Fourier-Taylor expansion 
# $$
# fₙ(θ) = ∑_{k=0}^{ n-1} f̂ₖ {\rm e}^{{\rm i}kθ}
# $$
# where $n = 20$,
# for $f(θ) = \exp(θ)$, $f(θ) = \exp(\cos(θ))$, and $\exp(\exp({\rm i} θ))$.
# For which of these functions does the approximation appear to converge?

n = 20
## TODO: plot the truncated Fourier-Taylor expansion with coefficients ranging from 0:n-1
## SOLUTION
f̂ = [fouriercoefficient(exp, k) for k = 0:n-1]
T = θ -> [exp(im*k*θ) for k = 0:n-1] # Create the Fourier basis
fₙ = θ -> transpose(T(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, exp.(g); label="f")
plot!(g, real(fₙ.(g)); label="Re fₙ(θ)")
plot!(g, imag(fₙ.(g)); label="Im fₙ(θ)")
## The approximation fails

f = θ -> exp(cos(θ))
f̂ = [fouriercoefficient(f, k) for k = 0:n-1]
fₙ = θ -> transpose(T(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, f.(g); label="f")
plot!(g, real(fₙ.(g)); label="Re fₙ(θ)")
plot!(g, imag(fₙ.(g)); label="Im fₙ(θ)")
## Fails to converge


f = θ -> exp(exp(im*θ))
f̂ = [fouriercoefficient(f, k) for k = 0:n-1]
fₙ = θ -> transpose(T(θ))*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, real.(f.(g)); label="Re f")
plot!(g, imag.(f.(g)); label="Im f")
plot!(g, real(fₙ.(g)); label="Re fₙ") 
plot!(g, imag(fₙ.(g)); label="Im fₙ")
## The real and imaginary parts match to high accuracy
## END

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
    ## SOLUTION
    ## We can either recall the integral formula or use the following to relate to the Fourier coefficient.
    ## We have
    ## f̂ₖ \exp(ikθ) = f̂ₖ \cos(kθ) + i f̂ₖ \sin(kθ) 
    ## f̂₋ₖ \exp(-ikθ) = f̂₋ₖ \cos(kθ) - i f̂₋ₖ \sin(kθ) 
    ## thus ignoring the sin terms we get f̌_k = f̂ₖ + f̂₋ₖ.
    ## apart from the `k = 0` case.
    ## Note in the case where f is even this simplifies as f̂ₖ = f̂₋ₖ.
    ## So your solution might be different if you used this property.
    if k ≠ 0
        real(fouriercoefficient(f, k) + fouriercoefficient(f, -k))
    else
        real(fouriercoefficient(f, 0))
    end
    ## END
end
@test cosinecoefficient(θ -> exp(cos(θ)), 5) isa Float64
@test cosinecoefficient(θ -> exp(cos(θ)), 5) ≈ 0.0005429263119137845

## TODO: plot the truncated cosine expansion with coefficients ranging from 0:n-1
## SOLUTION
f̂ = [cosinecoefficient(exp, k) for k = 0:n-1]
C = θ -> [cos(k*θ) for k = 0:n-1] # Create the Cosine basis
fₙ = θ -> C(θ)'f̂ # finite Cosine expansion, we can use ' as C(θ) is real now.
g = range(0, 2π, 1000) # plotting grid
plot(g, exp.(g); label="f")
plot!(g, fₙ.(g); label="fₙ(θ)")
## The approximation fails

f = θ -> exp(cos(θ))
f̂ = [cosinecoefficient(f, k) for k = 0:n-1]
fₙ = θ -> C(θ)'*f̂ # finite Fourier series
g = range(0, 2π, 1000) # plotting grid
plot(g, f.(g); label="f")
plot!(g, fₙ.(g); label="fₙ(θ)")
## Matches to high accuracy
## END


# ### V.1.1 Convergence of Fourier series

# Different function have different rates of convergence for their Fourier series,
# which is largely dictated by the rate of decay in their coefficients. We can explore
# how different _regularity_, i.e., the smoothness of the function, influences the rate of decay
# of a Fourier expansion.

# **Problem 2(a)** Plot  the absolute value of coefficients for the functions $\exp θ$, 
# $θ (2π - θ) \exp θ$ and $θ^2 (2π-θ)^2 \exp θ$. By scaling the $x$- and $y$-axis logarithmically deduce
# experimentally the rate of decay.

## TODO: plot the coefficients for the three functions with different smoothness properties and deduce the rate of decay.



# **Problem 2(b)** Repeat the above for $1/(\cos^2 θ + 1)$ and $1/(25 \cos^2 θ + 1)$, $\exp(\cos θ)$. Is the
# convergence algebraic, exponential, or super-exponential? 

## TODO: plot the coefficients for the three functions with different smoothness properties and deduce the rate of decay.



# ### V.1.2 Trapezium rule and discrete Fourier coefficients

# Above we use `quadgk` to compute the Fourier integrals, however, this has highly limited the number of coefficients
# we can compute as it becomes prohibitively expensive for large `k`. Instead, we can leverage the Trapezium rule which has
# highly structured error properties for computing Fourier integrals. We can use a modified version of the `trapeziumrule` function
# from Lab 1:

function periodictrapeziumrule(f, n)
    ret = 0.0
    for j = 0:n-1
        ret = ret + f(2π*j/n)
    end
    ret/n
end

f = θ -> exp(cos(θ))
@test periodictrapeziumrule(f, 10) ≈ quadgk(f, 0, 2π)[1]/(2π)

# **Problem 3** Implement the functio

# **Problem 4**


# ### V.1.3 Convergence of Approximate Fourier expansions


# **Problem 5**

# ## V.2 Discrete Fourier Transform


# The Fast Fourier Transform is 


# **Problem 6**

# ### V.2.1 Trigonometric Interpolation

# **Problem 7**