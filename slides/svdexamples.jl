using Images # Useful for plotting images given by pixels.
using MLDatasets: MNIST # load the MNIST function
mnist = MNIST() # loads the database
imgs, nums = mnist.features, mnist.targets

A = imgs[:,:,14004]

U,σ,V = svd(A)
r = 5

Gray.(A)
Gray.(U[:,1:r]*Diagonal(σ[1:r])*V[:,1:r]')

Gray.(U[:,1:r])
Gray.(Diagonal(σ[1:r]))
Gray.(V[:,1:r]')

length(U[:,1:r])+length(V[:,1:r]') # 280
length(A)