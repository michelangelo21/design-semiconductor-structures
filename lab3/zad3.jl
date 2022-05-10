using Plots
plotlyjs()

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

CsSiI₃ = @NamedTuple{Eg, Δ, γ₁, γ₂, γ₃, mₕ, Ep, a, α}((
    0.31, # Eg
    0.50, # Δ
    24.3, # γ₁
    11.5, # γ₂
    8.1, # γ₃
    unknown_mass_CsSiI3, # mₕ
    18.9, # Ep
    5.892, # a
    0.1, # α meV/K
))

CsPbₓSi₁₋ₓI₃(x, C=0) = @NamedTuple{Eg, Δ, γ₁, γ₂, γ₃, mₕ, Ep, a, α}(
    interpolation.(x, collect(CsPbI₃), collect(CsSiI₃), C)
)

x = 0.5
material = CsPbₓSi₁₋ₓI₃(x)

material.α
material.Eg

Eg(T, mat) = mat.Eg + mat.α * 1e-3 * T

VB₀ = 0.0
VB(T, mat) = VB₀
CS(T, mat) = VB₀ + Eg(T, mat)
CH(T, mat) = VB₀ + Eg(T, mat) + mat.Δ
CL(T, mat) = CH(T, mat)

Ts = 250:350 # K
plot(Ts, T -> VB(T, material), label="VB",
    xlabel="Temperature [K]",
    ylabel="Energy [eV]"
)
# plot(Ts, VB.(Ts, Ref(material)))
plot!(Ts, T -> CS(T, material), label="CS")
plot!(Ts, T -> CH(T, material), label="CH")
plot!(Ts, T -> CL(T, material), label="CL", linestyle=:dash)
