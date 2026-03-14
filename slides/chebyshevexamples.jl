using Plots, FFTW, LaTeXStrings

f = θ -> exp(θ) # function without a converging Fourier–Taylor expansion
g = range(-1,1,1000)
m = 10
n = 2m+1
θ = [2π/n*j for j=0:n-1]
𝐱 = θ/π .- 1
𝐟ₖ = fft(f.(𝐱))/n
fₙ = x -> transpose([exp(im*k*(π*(x+1))) for k=-m:m])*[𝐟ₖ[m+2:end]; 𝐟ₖ[1:m+1]]
plot(g, real.(f.(g)); label="exp(x)", linewidth=2, xlabel=L"x")
plot!(g, real.(fₙ.(g)); label=LaTeXString("\$f_{-$m:$m}(π(x+1))\$"), linewidth=2)
scatter!(𝐱, real.(f.(𝐱)); label=nothing) # we still interpolate exactly at the grid points
savefig("slides/exp.pdf")


f = θ -> exp(cos(θ)) # function without a converging Fourier–Taylor expansion
g = range(0,2π,1000)
m = 10
n = 2m+1
θ = [2π/n*j for j=0:n-1]
𝐟ₖ = fft(f.(θ))/n
fₙ = θ -> transpose([exp(im*k*(θ)) for k=-m:m])*[𝐟ₖ[m+2:end]; 𝐟ₖ[1:m+1]]
plot(g, real.(f.(g)); label="exp(cos(θ))", linewidth=2, xlabel=L"θ")
plot!(g, real.(fₙ.(g)); label=LaTeXString("\$f_{-$m:$m}(θ)\$"), linewidth=2)
scatter!(θ, real.(f.(θ)); label=nothing) # we still interpolate exactly at the grid points
savefig("slides/expcos.pdf")




f = x -> exp(x) # function without a converging Fourier–Taylor expansion
g = range(-1,1,1000)
m = 10
n = 2m+1
θ = [2π/n*j for j=0:n-1]
𝐟ₖ = fft(f.(cos.(θ)))/n
fₙ = θ -> transpose([exp(im*k*(θ)) for k=-m:m])*[𝐟ₖ[m+2:end]; 𝐟ₖ[1:m+1]]
plot(g, real.(f.(g)); label="exp(x)", linewidth=2, xlabel=L"x")
plot!(g, real.(fₙ.(acos.(g))); label="Chebyshev", linewidth=2)
scatter!(cos.(θ), real.(f.(cos.(θ))); label=nothing) # we still interpolate exactly at the grid points
savefig("slides/expcheb.pdf")



using ClassicalOrthogonalPolynomials

T = ChebyshevT()
# plot(T[:,1:5]; label=LaTeXString.("\$T_" .* string.((0:4)') .* "\$"), xlabel=L"x", linewidth=2)
plot(T[:,1:5]; legend=false, xlabel=L"x", linewidth=2)
savefig("slides/chebt.pdf")

g
plot(g, g .^(0:4)'; legend=false, xlabel=L"x", linewidth=2)
savefig("slides/monomial.pdf")


n = 10
plot(T[:,n+1]; label=L"T_{10}", xlabel=L"x", linewidth=2)
x = cos.([π/n*(j-1/2) for j=1:n])
scatter!(x, zero(x); label=nothing)
savefig("slides/chebtroots.pdf")