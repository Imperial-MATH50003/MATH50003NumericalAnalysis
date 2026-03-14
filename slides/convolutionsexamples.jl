using Images # Useful for plotting images given by pixels.
using MLDatasets: MNIST # load the MNIST function
mnist = MNIST() # loads the database
imgs, nums = mnist.features, mnist.targets

A = imgs[:,:,14004]
img = Gray.(A)
save("slides/default8.png", img)

m=3; n=3; img = Gray.(A[[end-m:end;1:end-m-1], [end-n:end;1:end-n-1]])
save("slides/shifted_default8.png", img)

m=11; n=14; img = Gray.(A[[m:end;1:m-1], [end-n:end;1:end-n-1]])
save("slides/periodic_default8.png", img)
