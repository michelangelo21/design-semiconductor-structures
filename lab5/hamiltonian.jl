# ? Jednostki Hartree ħ=1, m₀=1, e=1
using Plots
plotlyjs()

using LinearAlgebra
using Unitful, UnitfulAtomic
Unitful.preferunits(u"me_au, Å"...)
# promote energty to hartree by default
Unitful.promote_unit(::S, ::T) where {S<:Unitful.EnergyUnits,T<:Unitful.EnergyUnits} = u"Eh_au"



Base.@kwdef struct Perovskite
    Eg
    Δ
    γᴸ
    mₕ
    Ep
    a
    α
    c₁₁
    c₁₂
    aᵛ
    aᶜ
end


CsPbI₃ = Perovskite(;
    Eg=1.73u"eV", # experiment
    Δ=1.44u"eV",
    γᴸ=[9.1, 3.6, 0.7],
    mₕ=0.095u"me_au",
    Ep=41.6u"eV" / 2, # Ep = Eₚ_publikacja/2
    a=6.238u"Å",
    α=0.9u"meV/K",
    c₁₁=34.405u"GPa",
    c₁₂=4.709u"GPa",
    aᵛ=-2.762u"eV",
    aᶜ=-0.177u"eV"
)

CsSiI₃ = Perovskite(;
    Eg=0.31u"eV", # theory
    Δ=0.50u"eV",
    γᴸ=[24.3, 11.5, 8.1],
    mₕ=0.5 * (0.095 + 0.069) * u"me_au", # 0.5(CsPbI₃ + CsSnI₃)
    Ep=18.9u"eV" / 2, # Ep = Eₚ_publikacja/2
    a=5.892u"Å",
    α=0.7u"meV/K",
    # 1/3 ( CsPbI₃ + CsSnI₃ + CsGeI₃)
    c₁₁=(34.405 + 41.31 + 60.07) / 3 * u"GPa",
    c₁₂=(4.709 + 4.8 + 4.8) / 3 * u"GPa",
    aᵛ=(-2.762 + -3.651 + -2.257) / 3 * u"eV",
    aᶜ=(-0.177 + -0.052 + 0.971) / 3 * u"eV"
)

interpolation(x, q_A, q_B, C=zero(q_A)) = x * q_A + (1 - x) * q_B + x * (1 - x) * C

function interpolation(x, q_A::Perovskite, q_B::Perovskite) # TODO bowing
    Perovskite(
        (
            interpolation(x, getfield(q_A, f), getfield(q_B, f))
            for f in fieldnames(Perovskite)
        )...
    )
end

mat = interpolation(0.5, CsPbI₃, CsSiI₃)

"""
    band_structure(mat::Perovskite; T=300u"K", ϵₓₓ=0.0, h=0.01)

Compute band structure of perovskite on the path Γ -> R -> M.
knorm is used to scale the path.
knorm .< 0 means Γ -> R
knorm .> 0 means R -> M
# Examples
```julia-repl
julia> knorm, energies = band_structure(CsPbI₃; T=300u"K", ϵₓₓ=0.0)
"""
function band_structure(mat::Perovskite; T=300u"K", ϵₓₓ=0.0, knorm=-0.15:0.01:0.15)
    if !isa(T, Quantity)
        @error """Provide proper temperature units. e.g. T=300u"K" """
    end
    @assert isa(T, Quantity)

    Eg = mat.Eg + mat.α * T

    b = -2u"eV" # jak w GaAs

    ħ = 1u"ħ_au"
    m₀ = 1u"me_au"


    VBO = 0u"eV" |> auconvert # not VB0

    # bez naprężeń
    E_VB₀ = VBO
    E_CS₀ = E_VB₀ + Eg
    E_CH₀ = E_CS₀ + mat.Δ
    E_CL₀ = E_CH₀

    # naprężenia
    δE_Hⱽ = 2 * mat.aᵛ * (1 - mat.c₁₂ / mat.c₁₁) * ϵₓₓ
    δE_Hᶜ = 2 * mat.aᶜ * (1 - mat.c₁₂ / mat.c₁₁) * ϵₓₓ
    δE_S = b * (1 + mat.c₁₂ / mat.c₁₁) * ϵₓₓ

    # z naprężeniami
    VB₀ = E_VB₀ + δE_Hⱽ
    CS₀ = E_CS₀ + δE_Hᶜ
    CH₀ = E_CH₀ + δE_Hᶜ + δE_S
    CL₀ = E_CL₀ + δE_Hᶜ - δE_S




    # współczynniki z pub1
    γ = mat.γᴸ .- 1 ./ [3, 6, 6] * mat.Ep / Eg
    γᵥ = m₀ / mat.mₕ - mat.Ep / 3 * (2 / Eg + 1 / (Eg + mat.Δ))

    P = √(mat.Ep * ħ^2 / (2m₀)) |> auconvert


    # π/a * ( 0.15[1,1,1] -> [0,0,0] -> 0.15[0,0,1] )
    # knorm = -0.15:h:0.15
    ks = vcat(
        Ref([1, 1, 1]) .* knorm[knorm.<=0], # Γ - R
        Ref([0, 0, 1]) .* knorm[knorm.>0] # R - M
    ) .* (π / mat.a .|> auconvert)

    energies = map(ks) do k

        P₊ = P * (k[1] + im * k[2])
        P₋ = P * (k[1] - im * k[2])
        P_z = P * k[3]



        VB = VB₀ - ħ^2 / (2m₀) * γᵥ * (k[1]^2 + k[2]^2 + k[3]^2)
        CS = CS₀ + ħ^2 / (2m₀) * γ[1] * (k[1]^2 + k[2]^2 + k[3]^2)
        CH = CH₀ + ħ^2 / (2m₀) * ((γ[1] + γ[2]) * (k[1]^2 + k[2]^2) + (γ[1] - 2γ[2]) * k[3]^2)
        CL = CL₀ + ħ^2 / (2m₀) * ((γ[1] - γ[2]) * (k[1]^2 + k[2]^2) + (γ[1] + 2γ[2]) * k[3]^2)

        S = ħ^2 / (2m₀) * 2√3 * γ[3] * (-k[1] + im * k[2]) * k[3]
        R = ħ^2 / (2m₀) * √3 * (γ[2] * (k[1]^2 - k[2]^2) - 2im * γ[3] * k[1] * k[2])
        D = ħ^2 / (2m₀) * √2 * γ[2] * (k[1]^2 + k[2]^2 - 2k[3]^2)



        # TODO Hermitian from start
        # Hamiltonian
        H = [
            VB 0u"J" 1/√2*P₊ -√2/√3*P_z -1/√6*P₋ 0u"J" -1/√3*P_z -1/√3*P₋
            0u"J" VB 0u"J" 1/√6*P₊ -√2/√3*P_z -1/√2*P₋ -1/√3*P₊ 1/√3*P_z
            1/√2*P₋ 0u"J" CH S -R 0u"J" 1/√2*S -√2*R
            -√2/√3*P_z 1/√6*P₋ S' CL 0u"J" -R -D -√3/√2*S
            -1/√6*P₊ -√2/√3*P_z -R' 0u"J" CL -S -√3/√2*S' D
            0u"J" -1/√2*P₊ 0u"J" -R' -S' CH √2*R' 1/√2*S'
            -1/√3*P_z -1/√3*P₋ 1/√2*S' -D -√3/√2*S √2*R CS 0u"J"
            -1/√3*P₊ 1/√3*P_z -√2*R' -√3/√2*S' D 1/√2*S 0u"J" CS
        ]
        @assert H == H'
        Ĥ = Hermitian(H)

        Es = eigvals(Ĥ .|> austrip) * u"Eh_au" .|> u"eV"

        return Es
    end
    return knorm, energies
end

path_percent = 0.15
knorm, energies = band_structure(mat; knorm=-path_percent:0.005:path_percent)


function plot_energies(knorm, energies)
    plt = plot()
    for (i, lbl) in enumerate(reverse(["VB", "CS", "CH", "CL"]))
        plot!(plt, knorm, getindex.(energies ./ u"eV", 9 - 2i), label=lbl)
    end
    ylabel!(plt, "E [eV]")
    xlabel!(plt, "k_norm")
    path_percent = maximum(knorm)
    xticks!(plt, [-path_percent, 0, path_percent], ["$path_percent(Γ← R)", "R", "$path_percent(R→M)"])
    return plt
end

ϵₓₓs = -0.03:0.0002:0.03
@gif for ϵₓₓ in vcat(ϵₓₓs, reverse(ϵₓₓs))
    knorm, energies = band_structure(mat; ϵₓₓ=ϵₓₓ, knorm=-path_percent:0.01:path_percent)

    plt = plot_energies(knorm, energies)

    title!(plt, "ϵₓₓ = $(round(ϵₓₓ,digits=3))")
end

Ts = (250:1:350) * u"K"
@gif for T in vcat(Ts, reverse(Ts))
    knorm, energies = band_structure(mat; T=T, knorm=-path_percent:0.01:path_percent)

    plt = plot_energies(knorm, energies)

    title!(plt, "T = $T")
end

xs = 0:0.001:1
@gif for x in vcat(xs, reverse(xs))
    mat = interpolation(x, CsPbI₃, CsSiI₃)
    knorm, energies = band_structure(mat; knorm=-path_percent:0.01:path_percent)

    plt = plot_energies(knorm, energies)

    title!(plt, "x = $x")
end

[band_structure(mat; e)]
x = knorm
y = ϵₓₓs = -0.03:0.005:0.03 # ϵₓₓ

[]
surface(0:0.1:10, 0:0.1:10, foo)