using Plots

# linespacex = -8:0.01:8
# f(x) = (exp(x) + 1) / x
# plot(linespacex, f)
plot(-8:0.001:8, x -> (exp(x) + 1) / x)