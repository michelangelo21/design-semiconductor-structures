# Cs [Pb_x Si_{1-x}] I_3

using Plots
plotlyjs()

interpolation(x, q_A, q_B, C=0) = x * q_A + (1 - x) * q_B + x * (1 - x) * C

unknown_mass_CsSiI3 = 0.5 * (0.095 + 0.069)

CsPbI₃ = @NamedTuple{Eg, Δ, γ₁, γ₂, γ₃, mₕ, Ep, a}((
    1.73, # Eg
    1.44, # Δ
    9.1, # γ₁
    3.6, # γ₂
    0.7, # γ₃
    0.095, # mₕ
    41.6, # Ep
    6.238 # a
))

CsSiI₃ = @NamedTuple{Eg, Δ, γ₁, γ₂, γ₃, mₕ, Ep, a}((
    0.31, # Eg
    0.50, # Δ
    24.3, # γ₁
    11.5, # γ₂
    8.1, # γ₃
    unknown_mass_CsSiI3, # mₕ
    18.9, # Ep
    5.892 # a
))

CsPbₓSi₁₋ₓI₃(x, C=0) = @NamedTuple{Eg, Δ, γ₁, γ₂, γ₃, mₕ, Ep, a}(
    interpolation.(x, collect(CsPbI₃), collect(CsSiI₃), C)
)

xs = 0:0.01:1
ys = CsPbₓSi₁₋ₓI₃.(xs)
ys[1].Δ


plot(xs, [y.Eg for y in ys], xlabel=:x, ylabel="Eg [eV]", title=:CsPbₓSi₁₋ₓI₃)
savefig("interpolation_Eg.png")
plot(xs, [y.Δ for y in ys], xlabel=:x, ylabel="Δ [eV]", title=:CsPbₓSi₁₋ₓI₃)
savefig("interpolation_Δ.png")
plot(xs, [y.Ep for y in ys], xlabel=:x, ylabel="Ep [eV]", title=:CsPbₓSi₁₋ₓI₃)
savefig("interpolation_Ep.png")

plot(xs, [y.a for y in ys], xlabel=:x, ylabel="a [Å]", title=:CsPbₓSi₁₋ₓI₃)
savefig("interpolation_a.png")

for key in keys(ys[1])[3:6]
    plot(xs, [y[key] for y in ys], xlabel=:x, ylabel=key, title=:CsPbₓSi₁₋ₓI₃)
    savefig("interpolation_$key.png")
end

plot(xs, [y[gamma] for y in ys, gamma in [:γ₁, :γ₂, :γ₃]],
    xlabel=:x, ylabel=:γ, title=:CsPbₓSi₁₋ₓI₃)
savefig("interpolation_gammas.png")

plot([y.a for y in ys], [y.Eg for y in ys], xlabel="a [Å]", ylabel="Eg [eV]", title=:CsPbₓSi₁₋ₓI₃)
savefig("inter_a_Eg.png")


bowing = 1.5
ys_bow = CsPbₓSi₁₋ₓI₃.(xs, bowing)
plot(xs, [y.Eg for y in ys_bow], xlabel=:x, ylabel="Eg [eV]", title="CsPbₓSi₁₋ₓI₃, bowing = $bowing")
savefig("bowing15_Eg.png")
plot(xs, [y.Δ for y in ys_bow], xlabel=:x, ylabel="Δ [eV]", title="CsPbₓSi₁₋ₓI₃, bowing = $bowing")
savefig("bowing15_Δ.png")

plot([y.a for y in ys], [y.Eg for y in ys_bow], xlabel="a [Å]", ylabel="Eg [eV]", title="CsPbₓSi₁₋ₓI₃, bowing = $bowing")
savefig("bow15_a_Eg.png")