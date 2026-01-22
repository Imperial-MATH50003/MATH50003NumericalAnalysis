# # MATH50003 (2025–26)
# # Lab 2: I.3 Dual Numbers and I.4 Newton's Method

# In this lab we explore an alternative approach to computing derivatives:
# using _dual numbers_. This is a special mathematical object akin to complex numbers
# that allows us to compute derivatives to very high accuracy in an automated fashion.
# This is a basic example of [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)
# that is extremely important in Machine Learning and other computational applications.
# To realise dual numbers on a computer we need to introduce the notation of a _type_
# and create a customised type to represent dual numbers, which is what we discuss first.
# As an application of computing derivatives we consider root finding via [Newton's method](https://en.wikipedia.org/wiki/Newton%27s_method).


# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 1. Definition of dual numbers and functions applied dual numbers.
# 3. Newton's method for root finding.
#
# Coding knowledge:
#
# 1. The notion of a type and how to make your own type.
# 2. Defining functions whose arguments are restricted to specific types.
# 3. Overloading functions like `+`, `*`, and `exp` for a custom type.


# We load the `Test` and `Plots` packages to be used below.
# We also define the function `nanabs` which is useful for logarithmically scaled
# plots. For brevity we use a shorthand `cond ? expr1 : expr2` which just means
# ```julia
# if cond
#     expr1
# else
#     expr2
# end
# ```

using Plots, Test
nanabs(x) = x == 0 ? NaN : abs(x)


# ## Types in Julia


# Before we can use a concept like dual numbers we have to understand the notion of a _type_.
# In compiled languages like Julia everything has a type. The function `typeof` can be used to determine the type of,
# for example, a number.
# By default when we write an integer (e.g. `-123`) it is of type `Int`:

##

# On a 64-bit machine this will print `Int64`, where the `64` indicates it is using precisely 64 bits
# to represent the number (a topic we will come back to in Chapter II). If we write something with
# a decimal point it represents a "real" number, whose storage is of type `Float64`:

##

# This is called a _floating point number_, and again the `64` indicates it is using precisely
# 64 bits to represent this number. (We will see this is why computations like divided differences
# have large errors: because we are limiting the number of "digits" to represent numbers we need to
# round our computations.) Note that some operations involving `Int`s return `Float64`s:

##

# It is possible to have functions behave differently depending on the input type.
# To do so we can add a restriction denoted `::Int` or `::Float64` to the function _signature_.
# Here we create a function `foo` that is equal to `1` if the input is an `Int`, `0` if the input is
# a `Float64`, and `-1` otherwise:

##

# The last line returns a list of `Int`s, which has the type `Tuple`.
# Note that there is a difference between the set concept of integer and the type `Int`: whilst `3.0` is an integer
# its type is `Float64` so `foo(3.0) == 0`.

# **Remark** Every type has a _supertype_, which is an _abstract type_: something you can't make an instance of it.
# For example, in the same way that integers
# are subsets of the reals we have that `Int` and `Float64` are subtypes of
# `Real`. Which is a subtype of `Number`. Which, as is everything, a subtype of `Any`.

# Types allow for combining multiple numbers (or instances of other types) to represent a more complicated
# object. A simple example of this is a complex number,
# which stores two real numbers $x$ and $y$ (either `Int` or `Float64` or indeed other real number types not yet discussed)
# to represent the complex number $x + {\rm i} y$. In Julia ${\rm i} = \sqrt{-1}$ is denoted `im` and
# hence we can create a complex number like $1+2{\rm i}$ as follows:

##

# This complex number has two _fields_: the real and imaginary part. Accessing the fields is done
# using a `.`, here we display the real and imaginary parts as a `Tuple`:

##

# When we ask  its type we see it is a `Complex{Int}`:

##

# The `{Int}` indicates that each of the fields is an `Int`.
# Note we can add, subtract, multiply, or apply functions like `exp` to complex numbers:

##

# -----
# **Problem 1(a)** Use `typeof` to determine the type of `1.2 + 2.3im`.

## TODO: What is the type of 1.2 + 2.3im?



# **Problem 1(b)** Add another implementation of `foo` that returns `im` if the input
# is a `ComplexF64`.

## TODO: Overload foo for when the input is a ComplexF64 and return im


@test foo(1.1 + 2im) == im

# ------

# **Problem 2(a)** Consider the Taylor series approximation to the exponential:
# $$
# \exp z ≈ ∑_{k=0}^n {z^k \over k!}
# $$
# Complete the function `exp_t(z, n)` that computes this and returns a
# `Complex{Float64}` if the input is complex and a `Float64` if the input is real.
# Do not use the inbuilt `factorial` function.
# Hint: It might help to think inductively: for $s_k = z^k/k!$ we have
# $$
#   s_{k+1}  = {z \over k+1} s_k.
# $$

function exp_t(z, n)
    ## TODO: Compute the first (n+1)-terms of the Taylor series of exp
    ## evaluated at z
    
end

@test exp_t(1.0, 10) isa Float64 # isa is used to test the type of a result
@test exp_t(im, 10) isa ComplexF64 # isa is used to test the type of a result

@test exp_t(1.0, 100) ≈ exp(1)

# **Problem 2(b)** Plot the error for `n = 1:1000` of `exp_t(z, n)` for `z = 1, im, -5`, and `-100`,
# scaling the y-axis logarithmically.
# Does the method appear to converge for all values of $z$?

## TODO: plot the error for the Taylor series approximation.



# ------


# One of the powerful features of Julia is that it's very easy to make our own types. Let's begin with a simple
# implementation of a rational function $p/q$ where $p$ and $q$ are `Int`s.  Thus we want to create a new
# type called `Rat` with two fields `p` and `q` to represent the numerator and denominator, respectively.
# (For simplicity  we won't worry about restricting $p$ and $q$ to be `Int`.)
# We can construct such a type using the `struct` keyword:

##

# A new instance of `Rat` is created via `Rat(p,q)`, e.g., `Rat(1, 2)` represents `1/2`
# where the first argument specifies the numerator is `p = 1` and the second argument specifies the denominator is `q = 2`.
# The fields are accessed by `.`:

##

# Unfortunately we can't actually do anything with this type yet, so for example if we try to add two dual `Rat`s it will throw an error:

##

# The error is telling us to overload the `+` function when the inputs are both `Rat`.
# To do this we need to _import_ the `+` function and then we can overload it like any
# other function:

##

# We can support mixing `Rat` and `Int` by adding additional functionality:

##

# -----

# **Problem 3** Support `*`, `-`, `/`, and `==` for `Rat` and `Int`.

## We import +, -, *, / so we can _overload_ these operations
## specifically for Rat.
import Base: +, -, *, /, ==

## The ::Rat means the following version of == is only called if both
## arguments are Rat.
function ==(x::Rat, y::Rat)
    ## TODO: implement equality, making sure to check the case where
    ## the numerator/denominator are possibly reducible
    ## Hint: gcd and div may be useful. Use ? to find out what they do

    
end

## We can also support equality when x is a Rat and y is an Int
function ==(x::Rat, y::Int)
    ## TODO: implement
    
end

## TODO: implement ==(x::Int, y::Rat)


@test Rat(1, 2) == Rat(2, 4)
@test Rat(1, 2) ≠ Rat(1, 3)
@test Rat(2,2) == 1
@test 1 == Rat(2,2)

## TODO: implement +, -, *, and /,


@test Rat(1, 2) + Rat(1, 3) == Rat(5, 6)
@test Rat(1, 3) - Rat(1, 2) == Rat(-1, 6)
@test Rat(2, 3) * Rat(3, 4) == Rat(1, 2)
@test Rat(2, 3) / Rat(3, 4) == Rat(8, 9)

# ------

# ## I.3 Dual Numbers
#
# We now consider implementing a type `Dual` to represent the dual number $a + bϵ$,
# in a way similar to `Complex` or `Rat`. For simplicity we don't restrict the types of `a` and `b`
# but for us they will usually be `Float64`. We create this type very similar to `Rat` above:

##

# We can easily support addition of dual numbers as in `Rat` using the formula
# $$
# (a+bϵ) + (c+dϵ) = (a+c) + (b+d)ϵ
# $$

##

# For multiplication we used the fact that $ϵ^2 = 0$ to derive the formula
# $$
# (a+bϵ)*(c+dϵ) = ac +(bc+ad)ϵ.
# $$
# Here we support this operation by overloading `*` when the inputs are both
# `Dual`:

##


# ### I.3.1 Differentiating polynomials

# Dual numbers allow us to differentiate functions provided they are composed of
# operations overloaded for `Dual`. In particular, we have that
# $$
# f(x + b ϵ) = f(x) + bf'(x)ϵ
# $$
# and thus if we set `b = 1` the _dual part_ is equal to the derivative.
# We can use this fact to differentiate simple polynomials that only use `+`
# and `*`:

##

# A polynomial like `x^3 + 1` is not yet supported.
# To support this we need to add addition of `Dual` with `Int` or `Float64`.
# Note that both of these are _subtypes_ of `Real` and so restricting on `Real`
# will support both at the same time.
# We can overload the appropriate functions as follows:

##

# ### I.3.2 Differentiating functions

# We can also overload functions like `exp` so that they satisfy the rules of
# a _dual extension_, that is, are consistent with the formula $f(a+bϵ) = f(a) + bf'(a)ϵ$
# as follows:

##

# We can use this to differentiate a function that composes these basic operations:

##


# What makes dual numbers so effective is that, unlike divided differences, they are not
# prone to disasterous growth due to round-off errors: the above approximation
# matches the true answer to roughly 16 digits of accuracy.


# ------

# **Problem 4(a)** Add support for `-`, `cos`, `sin`, and `/` to the type `Dual`
# by replacing the `# TODO`s in the below code.


import Base: -, cos, sin, /

## The following supports negation -(a+bϵ)
-(x::Dual) = Dual(-x.a, -x.b)

## TODO: implement -(::Dual, ::Dual)



function cos(x::Dual)
    ## TODO: implement cos for Duals
    
end

function sin(x::Dual)
    ## TODO: implement sin for Duals
    
end

function /(x::Dual, y::Dual)
    ## TODO: implement division for Duals.
    ## Hint: think of this as x * (1/y)
    
end

x = 0.1
ϵ = Dual(0,1)
@test cos(sin(x+ϵ)/(x+ϵ)).b ≈ -((cos(x)/x - sin(x)/x^2)sin(sin(x)/x))


# **Problem 4(b)** Use dual numbers to compute the derivatives to
# 1. $\exp(\exp x \cos x + \sin x)$
# 2. $∏_{k=1}^{1000} \left({x \over k}-1\right)$
# 3. $f^{\rm s}_{1000}(x)$ where, as in Lab 1 Problem 3(d), $f^{\rm s}_n(x)$ corresponds to $n$-terms of the following continued fraction:
# $$
# 1 + {x-1 \over 2 + {x-1 \over 2 + {x-1 \over 2 + ⋱}}}.
# $$
# at the point $x = 0.1$. Compare with divided differences to give evidence that your implementation is correct.

## TODO: Use dual numbers to compute the derivatives of the 3 functions above.




# ------
# ## I.4 Newton's method

# Newton's method is a simple algorithmic approach that you may have seen before in school for computing roots (or zeros)
# of functions. The basic idea is given an initial guess $x_0$,
# find the first-order Taylor approximation $p(x)$ (i.e., find the line that matches the slope of the function at the point)
# $$
# f(x) ≈ \underbrace{f(x_0) + f'(x_0) (x- x_0)}_{p(x)}.
# $$
# We can then solve the root finding problem for $p(x)$ exactly:
# $$
# p(x) = 0 ⇔ x = x_0 - {f(x_0) \over f'(x_0)}.
# $$
# We take this root of $p(x)$ as the new initial guess and repeat. In other words, we have a simple sequence
# defined by
# $$
# x_{k+1} = x_k - {f(x_k) \over f'(x_k)}
# $$
# If the initial guess is "close enough" to a root $r$ of $f$ (ie $f(r) = 0$)
# then we have proven in the notes that $x_k → r$. Thus for large $N$ we have $x_N ≈ r$. 

# Dual numbers as implemented by `Dual` gives us a powerful tool to compute derivatives and get a simple implementation
# of Newton's method working:

##

# We can test that we have indeed found a root:
f(r)


# -----

# **Problem 5(a)** For $f(x) = x^5 + x^2 + 1$, plot the error of $x_k$ for `k = 1:15` where the
# y-axis is scaled logarithmically and chosen $x_0 = 0.1$ You may
# use the computed `r` as the "exact" root. What do you think the convergence rate is?

## TODO: compute and plot the error of `newton(f, 0.1, k)` for `k = 1:15`


# **Problem 5(b)** Use `newton` with a complex number to compute
# an approximation to a complex root of $f(x) = x^5 - x^2 + 1$.
# Verify the approximation is accurate by testing that it satisfies $f(r)$
# is approximately zero.


## TODO: By making the initial guess complex find a complex root.


# **Problem 5(c)** By changing the initial guesses compute 5 roots to
# $\sin x - 1/x$. Hint: you may need to add an overload for `/(x::Real, y::Dual)`.

## TODO: Use newton to compute roots of sin(x) - 1/x.

