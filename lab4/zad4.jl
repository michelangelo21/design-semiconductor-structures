using Plots
plotlyjs()

interpolation(x, q_A, q_B, C=0) = x * q_A + (1 - x) * q_B + x * (1 - x) * C

function interpolation(x, q_A::Dict, q_B::Dict, C=0)
    @assert keys(q_A) == keys(q_B)
    res = Dict()
    for (key, v_A) in q_A
        v_B = q_B[key]
        res[key] = interpolation(x, v_A, v_B, C)
    end
    return res
end

unknown_mass_CsSiI3 = 0.5 * (0.095 + 0.069)

CsPbI₃ = Dict(zip(["Eg", "Δ", "γ₁", "γ₂", "γ₃", "mₕ", "Ep", "a", "α"], [
    1.73, # Eg eV
    1.44, # Δ
    9.1, # γ₁
    3.6, # γ₂
    0.7, # γ₃
    0.095, # mₕ
    41.6, # Ep
    6.238, # a
    0.9, # α meV/K
]))

CsSiI₃ = Dict(zip(["Eg", "Δ", "γ₁", "γ₂", "γ₃", "mₕ", "Ep", "a", "α"], [
    0.31, # Eg
    0.50, # Δ
    24.3, # γ₁
    11.5, # γ₂
    8.1, # γ₃
    unknown_mass_CsSiI3, # mₕ
    18.9, # Ep
    5.892, # a
    0.1, # α meV/K
]))

# https://www.sciencedirect.com/science/article/pii/S003040261631021X?casa_token=o33T8v5jN3QAAAAA:I8IBYsFMKtk-aXxKz1GThGbgQfDKschBqnRRVHBQsx-ikcQH4iS-iVlYZ3DO8YZtA5JtWPR9eBo
CsPbI₃["c₁₁"] = 34.405 # GPa 
CsPbI₃["c₁₂"] = 4.709 # GPa

CsPbI₃["aᵛ"] = -2.762 # eV
CsPbI₃["aᶜ"] = -0.177 # eV

# CsSiI₃[val] = 1/3 ( CsPbI₃ + CsSnI₃ + CsGeI₃)
# https://www.worldscientific.com/doi/abs/10.1142/S0217984921500561 - CsSiI₃
# https://www.mdpi.com/2076-3417/10/15/5055 - CsGeI₃
CsSiI₃["c₁₁"] = (CsPbI₃["c₁₁"] + 41.31 + 60.07) / 3
CsSiI₃["c₁₂"] = (CsPbI₃["c₁₂"] + 3.69 + 48.61) / 3

CsSiI₃["aᵛ"] = (CsPbI₃["aᵛ"] + -3.651 + -2.257) / 3
CsSiI₃["aᶜ"] = (CsPbI₃["aᶜ"] + -0.052 + 0.971) / 3


function plot_energy_diagram_AB(x, percentage)
    b = -2 # like GaAs
    T = 300 # K

    mat_B = interpolation(x, CsPbI₃, CsSiI₃)
    mat_B["EgT"] = mat_B["Eg"] + mat_B["α"] * 1e-3 * T

    a_A = mat_B["a"] * (1 + percentage)
    ϵₓₓ = (a_A - mat_B["a"]) / mat_B["a"]

    δE_Hⱽ = 2 * mat_B["aᵛ"] * (1 - mat_B["c₁₂"] / mat_B["c₁₁"]) * ϵₓₓ
    δE_Hᶜ = 2 * mat_B["aᶜ"] * (1 - mat_B["c₁₂"] / mat_B["c₁₁"]) * ϵₓₓ
    δE_S = b * (1 + mat_B["c₁₂"] / mat_B["c₁₁"]) * ϵₓₓ

    VBO_A = 0 # eV
    VBO_B = 1 # eV
    Eg_A = VBO_B + mat_B["Eg"] + 3 # eV

    VB_B = VBO_B
    CS_B = VBO_B + mat_B["EgT"]
    CH_B = VBO_B + mat_B["EgT"] + mat_B["Δ"]
    CL_B = CH_B

    VB_A = VBO_A
    CS_A = VBO_A + Eg_A

    E_VB = VB_B .+ δE_Hⱽ
    E_CS = CS_B .+ δE_Hᶜ
    E_CH = CH_B .+ δE_Hᶜ .+ δE_S
    E_CL = CL_B .+ δE_Hᶜ .- δE_S

    d = 50 # nm
    X_A = vcat(-d:0, d:2d)
    X_B = 0:d

    fig = plot([-d:0; 0], [ones(d+1) * VB_A; E_VB], color=:blue, label="E_A_VB")
    plot!(fig, [d; d:2d], [E_VB; ones(d+1) * VB_A], color=:blue, label=:none)

    plot!(fig, [-d:0; 0], [ones(d + 1) * CS_A; E_CS], color=:orange, label="E_A_CS")
    plot!(fig, [d; d:2d], [E_CS; ones(d + 1) * CS_A], color=:orange, label=:none)

    plot!(fig, 0:d, ones(d + 1) * E_VB, label="E_B_VB", color=:blue, linestyle=:dash)
    plot!(fig, 0:d, ones(d + 1) * E_CS, label="E_B_CS", color=:orange, linestyle=:dash)
    plot!(fig, 0:d, ones(d + 1) * E_CH, label="E_B_CH", color=:red, linestyle=:dot)
    plot!(fig, 0:d, ones(d + 1) * E_CL, label="E_B_CL", color=:green, linestyle=:dash)

    title!(fig, "CsPbₓSi₁₋ₓI₃, x=$x, a_A=$(round(a_A,digits=3)), ϵₓₓ=$(round(ϵₓₓ*100,digits=3))%")
    xlabel!(fig, "[nm]")
    ylabel!(fig, "E [eV]")
    return fig
end
plot_energy_diagram_AB(0.8, 0.03)