# # MATH50003 (2023–24)
# # Lab 5: III.2 Differential Equations

# We also explore the reduction of differential equations to
# banded linear systems via divided differences. When we get lower bidiagonal systems these can be solved
# using forward substitution, whereas we will discuss the tridiagonal case later.

# We first load  packages we need including a couple new ones:


## LinearAlgebra contains routines for doing linear algebra
using LinearAlgebra, Plots, Test


# **Remark** One should normally not need to implement methods for solving differential equations
# oneself as there are packages available, including the high-performance
#  Julia package  [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Moreover Forward and Backward
# Euler are only the first baby steps to a wide range of time-steppers, with Runge–Kutta being
# one of the most successful.
# For example, in practice we can solve
# a simple differential equation like a pendulum $u'' = -\sin u$ can be solved
# as follows (writing at a system $u' = v, v' = -\sin u$):

using DifferentialEquations, LinearAlgebra, Plots

u = solve(ODEProblem((u,_,x) -> [u[2], -sin(u[1])], [1,0], (0,10)))
plot(u)

# However, even in these automated packages one has a choice of different methods with
# different behaviour, so it is important to understand on a mathematical level what is happening under the hood.


# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 2. Reduction of differential equations to bidiagonal or tridiagonal linear systems.
# 3. Two-point boundary value problems and their convergence rates.
#
# Coding knowledge:
#




#-----


# ## III.2 Differential Equations via Finite Differences

# We now turn to an important application of banded linear algebra:
# approximating solutions to linear differential equations. We will focus on first and second order
# but the techniques generalise beyond this, to vector problems, nonlinear differential equations, and partial differential equations.
# In particular we explore _finite difference_ approximations which use divided differences to replace derivatives.
# These are the most basic type of numerical method and many powerful alternatives
# exist, including Finite Element Methods and spectral methods.


# ### III.2.1 Indefinite integration

# We can use the right-sided divided difference to approximate derivatives.  Let's do an example of integrating $\cos x$ by discretising the ODE
# $$
#  u(0) = c, \qquad u'(x) = f(x)
# $$
# and see if our method matches
# the true answer of $\sin x$. Recall from the notes that this equation can be approximated by $u_k$ solving the bidiagonal linear system
# $$
# \begin{bmatrix}
#     1 \\
#     -1/h & 1/h \\
#     & ⋱ & ⋱ \\
#     && -1/h & 1/h \end{bmatrix} \begin{bmatrix}u_0\\u_1\\⋮\\u_n\end{bmatrix} = \begin{bmatrix}c\\ f(x_0)\\ ⋮ \\ f(x_{n-1})\end{bmatrix}.
# $$
# We can construct the bidiagonal matrix as follows:

n = 10
x = range(0, 1; length=n+1) # makes an n+1 point evenly spaced grid
h = step(x) # equivalent to 1/n
L = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)

# We can use this bidiagonal matrix along with `\` to solve the
# system via forward substitution:

c = 0 # u(0) = 0
f = x -> cos(x)

𝐟 = f.(x[1:end-1]) # evaluate f at all but the last point
𝐛 = [c; 𝐟]
𝐮 = L \ 𝐛 # integrate using forward-differences

plot(x, sin.(x); label="sin(x)", legend=:bottomright)
scatter!(x, 𝐮; label="forward")


#  We can estimate how fast it converges by measuring
# the ∞-norm error (using $\| 𝐱 \|_∞ := \max |x_k|$ which
# is implemented as `norm(x,Inf)`):

## Error from indefinite integration with c and f
function forward_err(u, c, f, n)
    x = range(0, 1; length = n+1)
    h = step(x) # equivalent to 1/n
    L = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
    𝐮 = L\ [c; f.(x[1:end-1])]
    errs = 𝐮 - u.(x) # compare numerics with "true" result
    norm(errs, Inf) # measure ∞-norm error
end


ns = 10 .^ (1:8) # solve up to n = 10 million
scatter(ns, forward_err.(sin, 0, f, ns); xscale=:log10, yscale=:log10, label="forward")
plot!(ns, ns .^ (-1); label="1/n", linestyle=:dash)

# We see that the method converges linearly (like $O(n^{-1})$).

# ------

# **Problem 5(a)** Implement Backward Euler as derived in the problem sheet to approximate
# indefinite-integration. How does the error compare to forward
# for $f(x) = \cos x$ on the interval $[0,1]$?
# Use the method to approximate the indefinite intergral of
# $$
# \exp(\exp x \cos x + \sin x)
# $$
# to 3 digits.

## TODO: Implement Backward Euler by constructing a lower bidiagonal linear system.
## SOLUTION

## To help with understanding we will also plot the errors even though this
## wasn't asked for in the question. We first compare backward and forward Euler:
c = 0 # u(0) = 0
f = x -> cos(x)
n = 10

x = range(0,1;length=n+1)
h=step(x)
A = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
ub = A\[c; f.(x[2:end])]
uf = A \ [c; f.(x[1:end-1])]

plot(x, sin.(x); label="sin(x)", legend=:bottomright)
scatter!(x, ub; label="backward")
scatter!(x, uf; label="forward")

#

## Comparing each method's errors, we see that the backward method has the same error as the forward method:


function forward_err(u, c, f, n)
    x = range(0, 1; length = n+1)
    h = step(x) # equivalent to 1/n
    L = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
    𝐮 = L\ [c; f.(x[1:end-1])]
    errs = 𝐮 - u.(x) # compare numerics with "true" result
    norm(errs, Inf) # measure ∞-norm error
end

function back_err(u, c, f, n)
    x = range(0,1;length=n+1)
    h=step(x)
    A = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
    ub = A\[c; f.(x[2:end])]
    norm(ub - u.(x), Inf)
end

c = 0 # u(0) = 0
f = x -> cos(x)
m = (x[1:end-1] + x[2:end])/2 # midpoints
ns = 10 .^ (1:8) # solve up to n = 10 million


scatter(ns, forward_err.(sin, 0, f, ns); xscale=:log10, yscale=:log10, label="forward")
scatter!(ns, back_err.(sin, 0, f, ns); label="back",alpha=0.5)
plot!(ns, ns .^ (-1); label="1/n")
plot!(ns, ns .^ (-2); label="1/n^2")


## Finally we get to the last part. One can increase n until
## 3 digits don't change. Or we can use the fact that the rate of convergence
## is O(n^(-1)) and guess that n = 100k suffices.


c = 0 # u(0) = 0
n = 100_000

##functions defined in the solutions to problem sheet 2
f = x -> exp(exp(x)cos(x) + sin(x))

x = range(0,1;length=n+1)
h=step(x)
A = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
uf = A\[c; f.(x[2:end])]
## Gives the values uf u on the grid to at least 3 digits

## END

# **Problem 5(b)** Implement indefinite-integration
# where we impose the equation on the midpoints $x̃_1,…,x̃_n$ defined as
# $$
# x̃_j = {x_{j+1} + x_j \over 2} = a + (j-1/2)h
# $$
# using the central difference formula
# $$
# u'(x̃_j) ≈ {u(x_j) - u(x_{j-1}) \over h}
# $$
# By plotting the errors show that this method converges at
# a faster rate than Forward or Backward Euler for $f(x) = \cos x$ on the interval $[0,1]$.


## TODO: Discretise at midpoints rather than our grid. The solution is still approximated on the original grid.
## SOLUTION


## The system is identical to before, just we evaluate the
## right-hand side at the midpoints.

n = 10
x = range(0, 1; length=n+1)
h = step(x)
A = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
c = 0 # u(0) = 0
f = x -> cos(x)

x̃ = (x[2:end] + x[1:end-1])/2 # could also be made with a comprehension
𝐟 = f.(x̃) # evaluate f at all but last points
𝐮 = A \ [c; 𝐟]

plot(x, sin.(x); label="sin(x)", legend=:bottomright)
scatter!(x, 𝐮; label="average grid point")


#


## Comparing the error to the midpoint method, we see that the errors are very similar:

function mid_err(u, c, f, n)
    x = range(0, 1; length = n+1)
    h = step(x) # equivalent to 1/n
    L = Bidiagonal([1; fill(1/h, n)], fill(-1/h, n), :L)
    x̃ = (x[2:end] + x[1:end-1])/2
    𝐟 = f.(x̃) # evaluate f at all but last points

    𝐮 = L\ [c; 𝐟]
    errs = 𝐮 - u.(x) # compare numerics with "true" result
    norm(errs, Inf) # measure ∞-norm error
end

c = 0 # u(0) = 0
f = x -> cos(x)
ns = 10 .^ (1:8) # solve up to n = 10 million


scatter(ns, mid_err.(sin, 0, f, ns); xscale=:log10, yscale=:log10, label="mid")
scatter!(ns, forward_err.(sin, 0, f, ns); label="forward")
plot!(ns, ns .^ (-1); label="1/n")
plot!(ns, ns .^ (-2); label="1/n^2")
## The error now decreases quadratically. Interestingly: we do not see
## the same growth in error as we did for computing derivatives. That is:
## solving an ODE is a more stable process than applying a differential operator.
## END

# ----

# ### III.2.2 Forward Euler


# We now adapt the approach for more general ODEs of the form
# $$
#   u'(x) + ω(x)u(x) = f(x), u(0) = c.
# $$
# We now have the system:
# $$
# \underbrace{\begin{bmatrix}
# 1 \\
# ω(x_0)-1/h & 1/h \\
# & ⋱ & ⋱ \\
# && ω(x_{n-1})-1/h & 1/h \end{bmatrix}}_L \underbrace{\begin{bmatrix}u_0 \\ u_1 \\ ⋮ \\ u_n\end{bmatrix} }_{𝐮} = \begin{bmatrix} c \\ f(x_0) \\ ⋮ \\ f(x_{n-1}) \end{bmatrix}
# $$
# Consider the simple example:
#  $$
#  u(0) = 1, u' + x u = {\rm e}^x
#  $$
#  which has an exact solution in terms of a special error function
#  (which I determined using Mathematica).


using SpecialFunctions
c = 1
ω = x -> x
n = 200
x = range(0, 1; length=n+1)
## exact solution, found in Mathematica
u = x -> -(1/2)*exp(-(1+x^2)/2)*(-2sqrt(ℯ) + sqrt(2π)erfi(1/sqrt(2)) - sqrt(2π)erfi((1 + x)/sqrt(2)))

h = step(x)
L = Bidiagonal([1; fill(1/h, n)], ω.(x[1:end-1]) .- 1/h, :L)

𝐛 = [c; exp.(x[1:end-1])]
𝐮 = L \ 𝐛

plot(x, u.(x); label="u(x)", legend=:bottomright)
scatter!(x, 𝐮; label="forward")

## We see that it is converging to the true result.

# ----


# **Problem  6** Implement backward Euler for solving:
# $$
# \begin{align*}
# u(0) &= 1, u'(t) - \cos(t) u(t) = t
# \end{align*}
# $$
# on the interval $[0,1]$. Approximate $u(1)$ to three digits accuracy.

## TODO: Implement backward Euler for the case with a variable coefficient.

## SOLUTION

function first_eq(n)
    x = range(0, 1; length=n+1)
    #find the step-size h
    h = step(x)
    L = Bidiagonal([1; fill(1/h, n) - cos.(x[2:end])], fill(-1/h, n), :L)
    L \ [1; x[2:end]]
end

for n in 2 .^ (1:13)
    println(first_eq(n)[end]) # print out final value
end
## We can guess that that $u(1) ≈ 2.96$ to three digits since those digits have stopped changing.

## END

# -----

# ### III.2.3 Poisson equation

# We now consider the Poisson equation with Dirichlet
# boundary conditions. In particular consider a case where
# we know the true answer: if $u(x) = \cos x^2$ then it solves the ODE:
# $$
# \begin{align*}
# u(0) = \underbrace{1}_c \\
# u''(x) = \underbrace{-4x^2 \cos(x^2) - 2\sin(x^2)}_{f(x)} \\
# u(1) = \underbrace{\cos 1}_d
# \end{align*}
# $$
# We approximate it by the solution to the tridiagonal system:
# $$
# \underbrace{\begin{bmatrix}
#     1 \\
#     1/h^2 & -2/h^2 & 1/h \\
#     & ⋱ & ⋱ & ⋱ \\
#    && 1/h^2 & -2/h^2 & 1/h \\
#    &&&& 1 \end{bmatrix}}_A \underbrace{\begin{bmatrix}u_0\\u_1\\⋮\\u_n\end{bmatrix} }_{𝐮} = \underbrace{\begin{bmatrix}c\\ f(x_0)\\ f(x_1)\\ ⋮ \\ f(x_{n-1})\\ d\end{bmatrix} }_{𝐛}
# $$
# We first construct the matrix $A$ using `Tridiagonal`:
n = 20
x = range(0, 1; length = n + 1)
h = step(x)
A = Tridiagonal([fill(1/h^2, n-1); 0],
                [1; fill(-2/h^2, n-1); 1],
                [0; fill(1/h^2, n-1)])

# Thus we get an approximation to our (known) solution:

u = x -> cos(x^2)
f = x -> -4x^2*cos(x^2) - 2sin(x^2)
𝐛 =  [1; f.(x[2:end-1]); cos(1)]
𝐮 = A \ 𝐛
plot(x, u.(x); label="u(x)", legend=:bottomright)
scatter!(x, 𝐮; label="finite differences")

# -----

# **Problem 7(a)** Estimate the rate of convergence in the ∞-norm using the previous example with an increasing number of grid points.

## TODO: Plot the ∞-norm error and estimate the convergence rate.
## SOLUTION

function poisson_err(u, c_0, c_1, f, n)
    x = range(0, 1; length = n+1)
    h = step(x)
    T = Tridiagonal([fill(1/h^2, n-1); 0], [1; fill(-2/h^2, n-1); 1], [0; fill(1/h^2, n-1)])
    uᶠ = T \ [c_0; f.(x[2:end-1]); c_1]
    norm(uᶠ - u.(x), Inf)
end



ns = 10 .^ (1:8) # solve up to n = 10 million
scatter(ns, poisson_err.(u, 1, cos(1), f, ns); xscale=:log10, yscale=:log10, label="error")
plot!(ns, ns .^ (-2); label="1/n^2")
## END


# **Problem 7(b)** Construct a finite-difference approximation to the
# forced Helmholtz equation
# $$
# \begin{align*}
# u(0) &= 0 \\
# u(1) &= 0 \\
# u'' + k^2 u &= {\rm e}^x
# \end{align*}
# $$
# and find an $n$ such  the error is less than $10^{-4}$ when compared
# with the true solution for $k=10$:
# $$
# u(x) = (-\cos(k x) + {\rm e}^x \cos(k x)^2 + \cot(k) \sin(k x) - {\rm e} \cos(k) \cot(k) \sin(k x) - {\rm e} \sin(k) \sin(k x) + {\rm e}^x \sin(k x)^2)/(1 + k^2)
# $$

## TODO: Generalise the second-order finite differences to allow for a $k^2 u$ term.

## SOLUTION
## We do something slightly different and use SymTridiagonal.
## You can also do this with Tridiagonal
function helm(k, n)
    x = range(0, 1; length = n+1)
    h = step(x)
    T = SymTridiagonal(ones(n-1)*(-2/h^2 + k^2),ones(n-2)*1/h^2)
    u = T \ exp.(x[2:end-1])
    [0; u; 0]
end

k = 10
u = x -> (-cos(k*x) + exp(x)cos(k*x)^2 + cot(k)sin(k*x) - ℯ*cos(k)cot(k)sin(k*x) - ℯ*sin(k)sin(k*x) + exp(x)sin(k*x)^2)/(1 + k^2)

n = 2048
x = range(0, 1; length=n+1)
@test norm(helm(k, n) - u.(x)) ≤ 1E-4
## END

