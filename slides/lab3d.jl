# # MATH50003 (2023–24)
# # Lab 3: II.1 Reals and II.2 Floating Point Arithmetic

# In this lab, we will explore how a computer represents integers (both signed and unsigned) and reals.
# In particular, its usage of modular and floating point arithmetic.
# We also investigate the rounding behaviour of floating point numbers
# and see how we can manually set the rounding. 



# **Learning Outcomes**
#
# Mathematical knowledge:
#
# 1. Behaviour of modular arithmetic for signed and unsigned integers.
# 2. Binary and hexadecimal integer formats.
# 3. Storage of floating point numbers.
#
# Coding knowledge:
#
# 1. Creating numbers with different formats via `UInt8`, `Float16`, etc.
# 2. The `sizeof`, `reinterpret`, `parse`, `typemax`, `bitstring`, and `printbits` functions
# 3. String construction and concatenation via `*`.
# 4. Creating floating point numbers by specifying their bits.
# 5. Setting the rounding mode in constructors like `Float32` and via `setrounding`.
# 6. High precision floating point numbers via `big` and setting precision via `setprecision`.


# We load an external `ColorBitstring` package
# which implements functions `printbits` (and `printlnbits`)
# to print the bits (and with a newline) of numbers in colour,
# and an external `SetRounding` package which allows for
# manually setting the rounding mode.

using ColorBitstring, SetRounding, Test

# If this fails you may need to call `] add ColorBitstring SetRounding`.

# ## Integers

# We first explore the representation and behaviour of integers in Julia,
# which is identical to other compiled languages like C, Swift, Rust, etc.
# and exposes how a computer actually computes with Integers.
# (Python is an interpreted language and its integers are more akin to
# `BigInt` discussed later, though in NumPy it uses the same
# types as Julia.)
# The basic integer type, which is constructed when we write say `5` is
# the signed integer type `Int`:

typeof(5)

# On a 64-bit machine this will return `Int64` indicating it is using 64-bits.
# (It's possible your machine is 32-bit in which case this will be `Int32`.)
# But other integer types exist, in particular for unsigned integers and for
# different number of bits.

# Every primitive number type is stored as a sequence of bits.
# The number of _bytes_ (i.e. 8-bits) can be deduced using the `sizeof` function:
##

# We will use the `printbits` command provided by `ColorBitstring` to see what bits
# are actually stored, eg.
##



# ### Unsigned integers

# Unsigned integers are used to represent non-negative integers. In Julia
# these correspond to types `UInt8`, `UInt16`, etc. where the number indicates
# the number of bits. The easiest way to create such an integer is to convert
# from an `Int`:

##

# This fails if a number cannot be represented as a specified type:
# e.g. `UInt8(-5)` and `UInt8(2^8)`.

# We can also create unsigned integers by specifying their bits
# by writing `0b` followed by a sequence of bits:

##
#
##


# -----
# **Problem 1(a)** Use binary format to create an `UInt32` corresponding to $(101101)_2$.

## TODO: Create an UInt32 representing (101101)_2



# **Problem 1(b)** What happens if you specify more than 64 bits using `0b⋅⋅…⋅⋅`?
# What if you specify more than 128 bits?

## TODO: Experiment with 0b with different amounts of digits.


# -----



# Integers use modular arithmetic for addition, subtraction and multiplication:

##



# If we go past the largest integer we overflow and wrap around:

##

# A similar phenomena happens with subtraction:

##


# Multiplication also works similarly. Multiplication by two shifts bits by
# one and modular arithmetic just drops extra bits so we have the following behaviour:

##

# ### Signed integers

# Signed integers represent negative and non-negative integers, where
# if the first bit is `1` it is interpreted as a negative number, according to
# the 2's complement format. There are multiple types of signed integers including
# `Int8`, `Int16`, `Int32`, and `Int64`. By default we create an `Int` but we can
# convert an `Int` to another signed integer type:
##

# It prints the same as `5` but calling `typeof` will confirm it is indeed an `Int8`.
# We can use `printbits` to see the expected binary format, matching that of `UInt8(5)`:
##

# Negative numbers use 2's complement. Suppose a
# $p$-bit signed integer has the same bits as a
# $p$-bit unsigned integer
# $x$. If The first bit is `0` this number is treated
# as positive and hence the signed integer is also $x$,
# i.e., `Int8(5)` and `Int8(5)` have the same bits.
# If the first bit of is `1` then it is interpreted as the
# negative number $x - 2^p$.
#
# For example, the bits `11111011` correspond to $251$
# but since the first bit is `1` it represents the integer
# $251 - 2^8 = -5$:
##

# We can use the `reinterpret` function to create an unsigned integer by
# specifying a sequence of bits and reinterpreting the bits as a signed integer:
##

# This is different from conversion via `Int8(0b11111111)` (which throws an error):
# `0b11111111` represents the (unsigned) integer $2^8-1 = 255$ and hence
# `Int8(0b11111111)` is equivalent to `Int8(UInt8(255))`. Since `255` is larger than
# the largest `Int8` (which is $2^7-1 = 127$) it would throw an error.


# -----

# **Problem 2** Construct an `Int8` representing the number $-63$ by specifying its
# bits and using `reinterpret`.

## TODO: Create an unsigned integer with 8 bits and reinterpret to get Int8(-63)

# -----

# ### Strings and parsing

# Strings are a convenient way of representing arbitrary strings of digits:
# we can convert bits of a number to a string of "1"s and "0"s
# using the function `bitstring`. For example:
##

# Whereas `printbits` prints the bits, this actually returns a string
# that can further be manipulated.


# We can `parse` a string of digits in base 2 or 10:
##

# Be careful with "negative" numbers, the following will fail: `parse(Int8, "10001011"; base=2)`

# It treats the string as binary digits, NOT bits. That is, negative numbers
# are represented using the minus sign:
##

# To concatenate strings we can use `*` (multiplication is used because string concatenation
# is non-commutative):
##

# The string consisting of the first nine characters can be found using `str[1:9]` where `str` is any string:
##

# The string consisting of the 11th through last character can be found using `str[11:end]`:
##

# We can combine string manipulation with `bitstring` and `parse` to manipulate bits.
# For example, we can see which `Int8` has the same bits as `Int8(5)` but with the third bit
# set to 1.
##


# -----

# **Problem 3(a)** Can you predict what the output of the following will be before hitting return?

bitstring(11)
#
bitstring(-11)



# **Problem 3(b)** Combine `parse`, `reinterpret`, and `UInt8` to convert the
# above string to a (negative) `Int8` with the specified bits.

## TODO: combine parse and reinterpret 


# **Problem 3(c)** Complete the following function that sets the 10th bit of an `Int32` to `1`,
# and returns an `Int32`, assuming that the input is a positive integer, using `bitstring`,
# `parse` and `*`.

function tenthbitto1(x::Int32)
    ## TODO: change the 10th bit of x to 1
    
end

@test tenthbitto1(Int32(100)) ≡ Int32(4194404)

# ------

# ### Hexadecimal and binary format

# In Julia unsigned integers are displayed in hexadecimal
# form: that is, in base-16.
# Since there are only 10 standard digits (`0-9`) it uses 6 letters (`a–f`) to represent
# 11–16. For example,
##
# because `f` corresponds to 15 and `a` corresponds to 10, and we have
# $$
# 15 * 16 + 10 = 250.
# $$
# The reason for this is that each hex-digit encodes 4 bits (since 4 bits have $2^4 = 16$ possible
# values) and hence two hex-digits are encode 1 byte, and thus the digits correspond
# exactly with how memory is divided into addresses.
# We can create unsigned integers either by specifying their hex format:
##
# Alternatively, we can specify their digits.
# For example, we know $(f)_{16} = 15 = (1111)_2$ and $(a)_{16} = 10 = (1010)_2$ and hence
# $250 = (fa)_{16} = (11111010)_2$ can be written as
##





# -----

# ## II.1 Reals
#
# Real numbers interpret a sequence of bits as a real number, specified in
# floating point.
# In Julia these correspond to 3 different floating-point types:

# • `Float64` is a type representing double precision ($F_{64} = F_{1023,11,52}$).
# We can create a `Float64` by including a
# decimal point when writing the number:

5.3 # isa Float64

# We can use `printbits` to see the stored bits:
##

# Or we can see the number in a more mathematical notation:
##

# The red bit is the sign bit (0 means positive). The Green bits represent the exponent, in this case:
##

# The blue bits are the significand. In this case represent $(1.01001)_2 = 1 + 2^{-2} + 2^{-5} = 1.28125$. And indeed
# we have
# $$
# 2^2 (1+2^{-2} + 2^{-5}) = 5 + 2^{-3} = 5.125
# $$

# Alternatively, one can use scientific notation: `5.125e0` to construct a `Float64`.
# `Float64` is the default format for
# scientific computing.
#
# • `Float32` is a type representing single precision ($F_{32} = F_{127,8,23}$).  We can create a `Float32` by including a
# `f0` when writing the number. Here we create a `Float32` and print its bits:
##

# Or more mathematically:
##

# Now the exponent is

σ = 127 # the shift according to Float32 format
0b10000001 - σ # == 2

# and again we see this represents `5.125`.
# `Float32` is generally the default format for graphics (on the _Graphics Processing Unit_, GPU),
# as the difference between 32 bits and 64 bits is indistinguishable to the eye in visualisation,
# and more data can be fit into a GPU's limited memory.

# • `Float16` is a type representing half-precision ($F_{16} = F_{15,5,10}$).
# It is important in machine learning where one wants to maximise the amount of data
# and high accuracy is not necessarily helpful.  We create `Float16` by converting a `Float64`
# as follows:
##

# I.e., mathematically:
##

# Now the exponent is

σ = 15 # the shift according to Float16 format
0b10001 - σ # == 2

# and we are representing $5.125$ again.
# `Float16` is important in machine learning where one wants to maximise the amount of data
# and high accuracy is not necessarily helpful.
#
# -----

# **Problem 4(a)** Use `reinterpret` and binary format to deduce which `Float16` has bits
# `1 01010 1010101010`.

## TODO: Construct a Float16 with the bits 1 01010 1010101010



#

# **Problem 4(b)** Use `printbits` to guess the binary representation of $1/5$.

## TODO: Use printbits and guess the binary expansion of 1/5.


# -----


# We now construct the largest and smallest `Float32` and check their bit sequences:
##

# For a given floating-point type, we can find these constants using the following functions:
##

# ### II.1.4 Sub-normal and special numbers

# If all the exponent bits are `0` then the number represents a "sub-normal" floating point number.

# **Example (creating a sub-normal number)** If we divide the smallest normal number by `2`, we get a subnormal number:
##

# Can you explain the bits?
#
# Zero is a sub-normal number, but it turns out there is also a negative zero:
##

#
# -----


# **Problem 5** Create the smallest positive non-zero sub-normal `Float16` by specifying
# its bits.

## TODO: create the smallest positive Float16


# -----

# The special numbers extend the real line by adding $±∞$ but also a notion of "not-a-number" ${\rm NaN}$.
# Whenever the bits of the exponent $q$ of a floating-point number are all 1 then they represent an element of $F^{\rm special}$.
# If all $b_k=0$, then the number represents either $±∞$, called `Inf` and `-Inf` for 64-bit floating-point numbers (or `Inf16`, `Inf32`
# for 16-bit and 32-bit, respectively):
##

# All other special floating-point numbers represent ${\rm NaN}$. One particular representation of ${\rm NaN}$
# is denoted by `NaN` for 64-bit floating-point numbers (or `NaN16`, `NaN32` for 16-bit and 32-bit, respectively):
##


# Arithmetic works differently on `Inf` and `NaN` and for undefined operations. 
# In particular we have:

1/0.0        # returns  Inf
1/(-0.0)     # returns -Inf
0.0/0.0      # returns  NaN
  
Inf*0        # returns  NaN
Inf+5        # returns  Inf
(-1)*Inf     # returns -Inf
1/Inf        # returns  0.0
1/(-Inf)     # returns -0.0
Inf - Inf    # returns  NaN
Inf ==  Inf  # returns  true
Inf == -Inf  # returns  false

NaN*0        # returns  NaN
NaN+5        # returns  NaN
1/NaN        # returns  NaN
NaN == NaN   # returns  false
NaN != NaN   # returns  true


# Essentially `NaN` is a CPU's way of indicating an error has occurred, but computation
# will continue.




# **Example (many `NaN`s)** What happens if we change some other $b_k$ to be nonzero?
# We can create bits as a string and see:

i = 0b0111110000010001 # an UInt16
reinterpret(Float16, i)

# Thus, there are many ways of representing `NaN`. (What a waste of perfectly good bit sequences!)

# -----

# ## II.2 Floating Point Arithmetic

# Real numbers that cannot be exactly represented by a specified floating point format are automatically rounded
# to the nearest float, however, the rouding mode can be changed to round up or down.
# In Julia, the rounding mode is specified by tags `RoundUp`, `RoundDown`, and
# `RoundNearest`. (There are also more exotic rounding strategies `RoundToZero`, `RoundNearestTiesAway` and
# `RoundNearestTiesUp` that we won't use.)



# Let's try rounding a `Float64` to a `Float32`.
##

# The default rounding mode can be changed:
##

# Or alternatively we can change the rounding mode for a chunk of code
# using `setrounding`. The following computes upper and lower bounds for `/`:
##


# **WARNING (compiled constants)**: Why did we first create a variable `x` instead of typing `1f0/3`?
# This is due to a very subtle issue where the compiler is _too clever for it's own good_: 
# it recognises `1f0/3` can be computed at compile time, but failed to recognise the rounding mode
# was changed. 

# **Problem 6** Complete functions `exp_t_3_down`/`exp_t_3_up` implementing the first
# three terms of the Taylor expansion of $\exp(x)$, that is, $1 + x + x^2/2 + x^3/6$ but where
# each operation is rounded down/up. Use `typeof(x)` to make sure you are changing the
# rounding mode for the right floating point type.

function exp_t_3_down(x)
    T = typeof(x) # use this to set the rounding mode
    ## TODO: use setrounding to compute 1 + x + x/2 + x^2/6 but rounding down
    
end

function exp_t_3_up(x)
    ## TODO: use setrounding to compute 1 + x + x/2 + x^2/6 but rounding up
    
end

@test exp_t_3_down(Float32(1)) ≡ 2.6666665f0 # ≡ checks type and all bits are equal
@test exp_t_3_up(Float32(1)) ≡ 2.6666667f0



# ### High-precision floating-point numbers


# It is possible to get higher precision (more signficand and exponent bits)
#  of a floating-point number
# using the `BigFloat` type, which results from the usage of `big`
# when the result is not an integer.
# For example, here is an approximation of 1/3 accurate
# to 77 decimal digits:
##

# We can see a mathematical version of what's stored:
##

# Note we can set the rounding mode as in `Float64`, e.g., 
# this gives (rigorous) bounds on
# `1/3`:
##

# We can also increase the precision, e.g., this finds bounds on `1/3` accurate to 
# more than 1000 decimal places:
##


# **Problem 7** Inbuilt functions like `exp`, `sqrt`, etc. support `BigFloat`.
# Compute at least the first thousand decimal digits of `ℯ` using `setprecision`
# and the inbuilt `exp` function.

## TODO: Use big and setprecision to compute the first thousand digits of ℯ.


