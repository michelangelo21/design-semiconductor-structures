interpolation(x, q_A, q_B, C=0) = x * q_A + (1 - x) * q_B + x * (1 - x) * C

unknown_mass_CsSiI3 = 0.5 * (0.095 + 0.069)

CsPbI₃ = @NamedTuple{Eg, Δ, γ₁, γ₂, γ₃, mₕ, Ep, a, α}((
    1.73, # Eg eV
    1.44, # Δ
    9.1, # γ₁
    3.6, # γ₂
    0.7, # γ₃
    0.095, # mₕ
    41.6, # Ep
    6.238, # a
    0.9, # α meV/K
))

c₁₁ = 1
c₁₂ = 2 * c₁₁
aⱽ_Pb = -2.762
aᶜ_Pb = -0.177

aⱽ = aⱽ_Pb
aᶜ = aᶜ_Pb
b = -2
ϵₓₓ = -0.02:0.001:0.02

δE_Hⱽ = 2 * aⱽ * (1 - c₁₂ / c₁₁) * ϵₓₓ
δE_Hᶜ = 2 * aᶜ * (1 - c₁₂ / c₁₁) * ϵₓₓ
δE_S = b * (1 + c₁₂ / c₁₁) * ϵₓₓ

Eg = CsPbI₃.Eg
Δ = CsPbI₃.Δ

VB₀ = 0.0
VB = VB₀
CS = VB₀ + Eg
CH = VB₀ + Eg + Δ
CL = CH

E_VB = VB .+ δE_Hⱽ
E_CS = CS .+ δE_Hᶜ
E_CH = CH .+ δE_Hᶜ .+ δE_S
E_CL = CL .+ δE_Hᶜ .- δE_S

using Plots
plotlyjs()
plot(ϵₓₓ, E_VB, label="E_VB")
plot!(ϵₓₓ, E_CS, label="E_CS")
plot!(ϵₓₓ, E_CH, label="E_CH")
plot!(ϵₓₓ, E_CL, label="E_CL")

