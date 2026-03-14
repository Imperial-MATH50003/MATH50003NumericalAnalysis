using Plots, FFTW

f = θ -> exp(θ) # function without a converging Fourier–Taylor expansion
g = range(0,2π,1000)
plot(g, real.(f.(g)); label="exp(θ)")
n = 4
θ = [2π/n*j for j=0:n-1]
𝐟ₖ = fft(f.(θ))/n
fₙ = θ -> transpose([exp(im*k*θ) for k=0:n-1])𝐟ₖ
plot!(g, real.(fₙ.(g)); label="n = $n")
scatter!(θ, real.(f.(θ)); label=nothing) # we still interpolate exactly at the grid points