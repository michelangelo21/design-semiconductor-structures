using Plots
plotlyjs()

interpolation(x, q_A, q_B, C=0) = x * q_A + (1 - x) * q_B + x * (1 - x) * C

# Ga_x In_(1-x) | P
Eg_GaP = 2.886 # eV
Eg_InP = 1.4236 # eV
C = 0.65 # bowing

Eg_GaInP(x) = interpolation(x, Eg_GaP, Eg_InP, C)

plot(0:0.01:1, Eg_GaInP, xlabel="x", ylabel="E_g [eV]", title="Gaₓ In₍₁₋ₓ₎ P")
