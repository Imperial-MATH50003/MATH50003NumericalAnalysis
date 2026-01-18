# # MATH50003 (2024–25)
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
## DEMO
sizeof(Int64) # 8 bytes == 8*8 bits == 64 bits
## END

# We will use the `printbits` command provided by `ColorBitstring` to see what bits
# are actually stored, eg.
## DEMO
printbits(5)
## END



# ### Unsigned integers

# Unsigned integers are used to represent non-negative integers. In Julia
# these correspond to types `UInt8`, `UInt16`, etc. where the number indicates
# the number of bits. The easiest way to create such an integer is to convert
# from an `Int`:

## DEMO
UInt8(5) # creates an Int and converts it to an UInt8
         ## displaying the result in hex, i.e. base-16
## END

# This fails if a number cannot be represented as a specified type:
# e.g. `UInt8(-5)` and `UInt8(2^8)`.

# We can also create unsigned integers by specifying their bits
# by writing `0b` followed by a sequence of bits:

## DEMO
0b101 # creates an UInt8, the smallest type with at least 3 bits
## END
#
## DEMO
0b10111011101 # creates an UInt16, the smallest type with at least 11 bits
## END


# -----
# **Problem 1(a)** Use binary format to create an `UInt32` corresponding to $(101101)_2$.

## TODO: Create an UInt32 representing (101101)_2
## SOLUTION
UInt32(0b101101) # without the UInt32 it will be a UInt8. Another solution would be 0b00000000000101101
## END


# **Problem 1(b)** What happens if you specify more than 64 bits using `0b⋅⋅…⋅⋅`?
# What if you specify more than 128 bits?

## TODO: Experiment with 0b with different amounts of digits.
## SOLUTION
typeof(0b111111111111111111111111111111111111111111111111111111111111111111111111111111111111) # creates a UInt128

typeof(0b111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111) # creates a BigInt
## END

# -----



# Integers use modular arithmetic for addition, subtraction and multiplication:

## DEMO
x = UInt8(17)  # An 8-bit representation of the number 17, i.e. with bits 00010001
y = UInt8(3)   # An 8-bit representation of the number   3, i.e. with bits 00000011
printbits(x); println(" + "); printbits(y); println(" = ")
printlnbits(x + y) # + is automatically modular arithmetic
printbits(x); println(" - "); printbits(y); println(" = ")
printbits(x - y) # - is automatically modular arithmetic
## END



# If we go past the largest integer we overflow and wrap around:

## DEMO
x = UInt8(255) # An 8-bit representation of the number 255, i.e. with bits 11111111
y = UInt8(1)   # An 8-bit representation of the number   1, i.e. with bits 00000001
printbits(x); println(" + "); printbits(y); println(" = ")
printbits(x + y) # + is automatically modular arithmetic
## END

# A similar phenomena happens with subtraction:

## DEMO
x = UInt8(3) # An 8-bit representation of the number   3, i.e. with bits 00000011
y = UInt8(5) # An 8-bit representation of the number   5, i.e. with bits 00000101
printbits(x); println(" - "); printbits(y); println(" = ")
printbits(x - y) # - is automatically modular arithmetic
## END


# Multiplication also works similarly. Multiplication by two shifts bits by
# one and modular arithmetic just drops extra bits so we have the following behaviour:

## DEMO
x = UInt8(254) # An 8-bit representation of the number 254, i.e. with bits 11111110
y = UInt8(2)   # An 8-bit representation of the number   2, i.e. with bits 00000010
printbits(x); println(" * "); printbits(y); println(" = ")
printbits(x * y) # represents 252
## END

# ### Signed integers

# Signed integers represent negative and non-negative integers, where
# if the first bit is `1` it is interpreted as a negative number, according to
# the 2's complement format. There are multiple types of signed integers including
# `Int8`, `Int16`, `Int32`, and `Int64`. By default we create an `Int` but we can
# convert an `Int` to another signed integer type:
## DEMO
Int8(5) # display of Int8 does not reveal its type
## END

# It prints the same as `5` but calling `typeof` will confirm it is indeed an `Int8`.
# We can use `printbits` to see the expected binary format, matching that of `UInt8(5)`:
## DEMO
printbits(Int8(5)) # 5 = 2^2 + 1
## END

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
## DEMO
printbits(Int8(-5))
## END

# We can use the `reinterpret` function to create an unsigned integer by
# specifying a sequence of bits and reinterpreting the bits as a signed integer:
## DEMO
reinterpret(Int8, 0b11111111) # Create an Int8 with the bits 11111111
## END

# This is different from conversion via `Int8(0b11111111)` (which throws an error):
# `0b11111111` represents the (unsigned) integer $2^8-1 = 255$ and hence
# `Int8(0b11111111)` is equivalent to `Int8(UInt8(255))`. Since `255` is larger than
# the largest `Int8` (which is $2^7-1 = 127$) it would throw an error.


# -----

# **Problem 2** Construct an `Int8` representing the number $-63$ by specifying its
# bits and using `reinterpret`.

## TODO: Create an unsigned integer with 8 bits and reinterpret to get Int8(-63)
## SOLUTION
## -63 + 256 = 193 = 128 + 64 + 1 = 2^7 + 2^6 + 1 = (11000001)_2

reinterpret(Int8,0b11000001)

## END
# -----

# ### Strings and parsing

# Strings are a convenient way of representing arbitrary strings of digits:
# we can convert bits of a number to a string of "1"s and "0"s
# using the function `bitstring`. For example:
## DEMO
bitstring(Int8(5))
## END

# Whereas `printbits` prints the bits, this actually returns a string
# that can further be manipulated.


# We can `parse` a string of digits in base 2 or 10:
## DEMO
parse(Int8, "11"; base=2), # represents 2 + 1 = 3 as an Int8
parse(Int8, "00001011"; base=2) # represents 2^3 + 2 + 1 = 11 as an Int8
## END

# Be careful with "negative" numbers, the following will fail: `parse(Int8, "10001011"; base=2)`

# It treats the string as binary digits, NOT bits. That is, negative numbers
# are represented using the minus sign:
## DEMO
parse(Int8, "-00001011"; base=2)
## END

# To concatenate strings we can use `*` (multiplication is used because string concatenation
# is non-commutative):
## DEMO
"hi" * "bye"
## END

# The string consisting of the first nine characters can be found using `str[1:9]` where `str` is any string:
## DEMO
str = "hibye0123445556"
str[1:9]  # returns "hibye0123"
## END

# The string consisting of the 11th through last character can be found using `str[11:end]`:
## DEMO
str = "hibye0123445556"
str[11:end]  # returns "45556"
## END

# We can combine string manipulation with `bitstring` and `parse` to manipulate bits.
# For example, we can see which `Int8` has the same bits as `Int8(5)` but with the third bit
# set to 1.
## DEMO
str = bitstring(Int8(5)) # string of bits for 5, eg  "00000101"
tru = str[4:end] # drop first four characters of the string, eg "000101"
swa = str[1:2] * "1" * tru # add the character "1" at the third position, eg "00100101"
parse(Int8, swa; base=2) # answer is 37 = 5 + 2^5
## END


# -----

# **Problem 3(a)** Can you predict what the output of the following will be before hitting return?

bitstring(11)
#
bitstring(-11)

## SOLUTION
bitstring(11) # "0000000000000000000000000000000000000000000000000000000000001011"
bitstring(-11) # "1111111111111111111111111111111111111111111111111111111111110101"
## this is because mod(-11, 2^64) == 2^64 - 12 == 0b10000…000 - 0b1100 == 0b111…11 - 0b1011 + 0b1
## END

# **Problem 3(b)** Combine `parse`, `reinterpret`, and `UInt8` to convert the
#  string  `"10001011"` to a (negative) `Int8` with the specified bits.

## TODO: combine parse and reinterpret 
## SOLUTION

## The  code `parse(Int8, "-00001011"; base=2)` creates an `Int8` with bits "11110101". Instead, we first parse the bits:

x = reinterpret(Int8, parse(UInt8, "10001011"; base=2)) # -117
bitstring(x)

## END

# **Problem 3(c)** Complete the following function that sets the 10th bit of an `Int32` to `1`,
# and returns an `Int32`, assuming that the input is a positive integer, using `bitstring`,
# `parse` and `*`.

function tenthbitto1(x::Int32)
    ## TODO: change the 10th bit of x to 1
    ## SOLUTION
    ret = bitstring(x)
    parse(Int32, ret[1:9] * "1" * ret[11:end]; base=2)
    ## END
end

@test tenthbitto1(Int32(100)) ≡ Int32(4194404)

# ------

# ### Hexadecimal and binary format

# In Julia unsigned integers are displayed in hexadecimal
# form: that is, in base-16.
# Since there are only 10 standard digits (`0-9`) it uses 6 letters (`a–f`) to represent
# 11–16. For example,
## DEMO
UInt8(250)
## END
# because `f` corresponds to 15 and `a` corresponds to 10, and we have
# $$
# 15 * 16 + 10 = 250.
# $$
# The reason for this is that each hex-digit encodes 4 bits (since 4 bits have $2^4 = 16$ possible
# values) and hence two hex-digits are encode 1 byte, and thus the digits correspond
# exactly with how memory is divided into addresses.
# We can create unsigned integers either by specifying their hex format:
## DEMO
0xfa
## END
# Alternatively, we can specify their digits.
# For example, we know $(f)_{16} = 15 = (1111)_2$ and $(a)_{16} = 10 = (1010)_2$ and hence
# $250 = (fa)_{16} = (11111010)_2$ can be written as
## DEMO
0b11111010
## END





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
## DEMO
printbits(5.125)
## END

# Or we can see the number in a more mathematical notation:
## DEMO
binarystring(5.125)
## END

# The red bit is the sign bit (0 means positive). The Green bits represent the exponent, in this case:
## DEMO
σ = 1023 # the shift according to Float64 format
0b10000000001 - σ # == 2
## END

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
## DEMO
printbits(5.125f0) # 5.125f0 of type Float32
## END

# Or more mathematically:
## DEMO
binarystring(5.125f0)
## END

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
## DEMO
printbits(Float16(5.125))
## END

# I.e., mathematically:
## DEMO
binarystring(Float16(5.125))
## END

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
## SOLUTION
reinterpret(Float16, 0b1010101010101010)
## END


#

# **Problem 4(b)** Use `printbits` to guess the binary representation of $1/5$.

## TODO: Use printbits and guess the binary expansion of 1/5.
## SOLUTION

printbits(1/5)
## exponent is 0b01111111100 == 1020 so we have 2^(1020 - 1023) = 2^(-3)
## significand is 1.1001100110011001100110011001100110011001100110011010
## guess: 1/5 == 2^(-3) (1.10011001100…)_2 2^(-3) (∑_{k=0}^∞ (2^(-4k) + 2^(-4k-1)))

## END

# -----


# We now construct the largest and smallest `Float32` and check their bit sequences:
## DEMO
σ,Q,S = 127,8,23 # Float32
εₘ = 2.0^(-S)
printlnbits(Float32(2.0^(1-σ))) # smallest positive normal Float32
printlnbits(Float32(2.0^(2^Q-2-σ) * (2-εₘ))) # largest normal Float32
## END

# For a given floating-point type, we can find these constants using the following functions:
## DEMO
eps(Float32), floatmin(Float32), floatmax(Float32)
## END

# ### II.1.4 Sub-normal and special numbers

# If all the exponent bits are `0` then the number represents a "sub-normal" floating point number.

# **Example (creating a sub-normal number)** If we divide the smallest normal number by `2`, we get a subnormal number:
## DEMO
mn = floatmin(Float32) # smallest normal Float32
printlnbits(mn)
printbits(mn/2)
## END

# Can you explain the bits?
#
# Zero is a sub-normal number, but it turns out there is also a negative zero:
## DEMO
printlnbits(0.0) # 0 has all bits 0
printlnbits(-0.0) # -0 has sign bit 1 and all other bits zero
## END

#
# -----


# **Problem 5** Create the smallest positive non-zero sub-normal `Float16` by specifying
# its bits.

## TODO: create the smallest positive Float16
## SOLUTION
## sign is + so sign bit is 0, exponent is 00000 and significand is all zeros apart from a 1:
reinterpret(Float16, 0b0000000000000001) # == nextfloat(Float16(0))
## END

# -----

# The special numbers extend the real line by adding $±∞$ but also a notion of "not-a-number" ${\rm NaN}$.
# Whenever the bits of the exponent $q$ of a floating-point number are all 1 then they represent an element of $F^{\rm special}$.
# If all $b_k=0$, then the number represents either $±∞$, called `Inf` and `-Inf` for 64-bit floating-point numbers (or `Inf16`, `Inf32`
# for 16-bit and 32-bit, respectively):
## DEMO
printlnbits(Inf16)
printbits(-Inf16)
## END

# All other special floating-point numbers represent ${\rm NaN}$. One particular representation of ${\rm NaN}$
# is denoted by `NaN` for 64-bit floating-point numbers (or `NaN16`, `NaN32` for 16-bit and 32-bit, respectively):
## DEMO
printbits(NaN16)
## END


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
## DEMO
printlnbits(1/3)  # 64 bits
printbits(Float32(1/3))  # round to nearest 32-bit
## END

# The default rounding mode can be changed:
## DEMO
printbits(Float32(1/3,RoundDown)) # Rounds from a Float64 to Float32, rounding down
## END

# Or alternatively we can change the rounding mode for a chunk of code
# using `setrounding`. The following computes upper and lower bounds for `/`:
## DEMO
x = 1f0
setrounding(Float32, RoundDown) do
  x/3
end,
setrounding(Float32, RoundUp) do
  x/3
end
## END


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
    ## SOLUTION
    setrounding(T, RoundDown) do
        1 + x + x^2/2 + x^3/6
    end
    ## END
end

function exp_t_3_up(x)
    ## TODO: use setrounding to compute 1 + x + x/2 + x^2/6 but rounding up
    ## SOLUTION
    T = typeof(x) # use this to set the rounding mode
    setrounding(T, RoundUp) do
        1 + x + x^2/2 + x^3/6
    end
    ## END
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
## DEMO
big(1)/3
## END

# We can see a mathematical version of what's stored:
## DEMO
binarystring(big(1)/3)
## END

# Note we can set the rounding mode as in `Float64`, e.g., 
# this gives (rigorous) bounds on
# `1/3`:
## DEMO
setrounding(BigFloat, RoundDown) do
  big(1)/3
end, setrounding(BigFloat, RoundUp) do
  big(1)/3
end
## END

# We can also increase the precision, e.g., this finds bounds on `1/3` accurate to 
# more than 1000 decimal places:
## DEMO
setprecision(4_000) do # 4000 bit precision
  setrounding(BigFloat, RoundDown) do
    big(1)/3
  end, setrounding(BigFloat, RoundUp) do
    big(1)/3
  end
end
## END


# **Problem 7** Inbuilt functions like `exp`, `sqrt`, etc. support `BigFloat`.
# Compute at least the first thousand decimal digits of `ℯ` using `setprecision`
# and the inbuilt `exp` function.

## TODO: Use big and setprecision to compute the first thousand digits of ℯ.
## SOLUTION
x = setprecision(4_000) do
    exp(big(1.0))
end

length(string(x)) == 1207 # we have 1205 digits

## END

