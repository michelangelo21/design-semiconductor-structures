using Plots
plotlyjs()

q_AB(x, q_A, q_B, C=0) = x * q_A + (1 - x) * q_B + x * (1 - x) * C

# Ga_x In_(1-x) | P
Eg_GaP = 2.886 # eV
Eg_InP = 1.4236 # eV
C = 0.65

Eg_GaInP(x) = q_AB(x, Eg_GaP, Eg_InP, C)

plot(0:0.01:1, Eg_GaInP)
