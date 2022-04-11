using Plots

import PlotlyJS

plotlyjs()

unknown_mass_CsSiI₃ = 0.5 * (0.095 + 0.069) # średnia z CsXI₃

struct Perovskite{T}
    Eg::T
    Δ::T
    γ₁::T
    γ₂::T
    γ₃::T
    mₕ::T
    Ep::T
    a::T
end

CsPbI₃ = Perovskite(
    1.73, # Eg
    1.44, # Δ
    9.1, # γ₁
    3.6, # γ₂
    0.7, # γ₃
    0.095, # mₕ
    41.6, # Ep
    6.238 # a
)

CsSiI₃ = Perovskite(
    0.31, # Eg
    0.50, # Δ
    24.3, # γ₁
    11.5, # γ₂
    8.1, # γ₃
    unknown_mass_CsSiI₃, # mₕ
    18.9, # Ep
    5.892 # a
)

interpolation(x, q_A, q_B; C=0) = x * q_A + (1 - x) * q_B + x * (1 - x) * C


function interpolation(x, pA::Perovskite, pB::Perovskite, C=0)

end

tmp = CsPbI₃
propertynames(1.0)

# ╔═╡ 831bf8a4-0e9a-4fd8-87ac-8c1ce88db102
CsPbₓSi₁₋ₓI₃(x; C=0) = Perovskite(
    interpolation.(x, CsPbI₃, CsSiI₃; C=C)
)
interpolation.(0.7, CsPbI₃, CsSiI₃; C=0)
CsPbI₃ .+ CsSiI₃
Broadcast.broadcastable(CsPbI₃)
Ref(CsPbI₃)
Broadcast.broadcastable(p::Perovskite) = Ref(p)
collect(CsPbI₃)
@benchmark Perovskite((getfield($tmp, prop) for prop in propertynames($CsPbI₃))...)
fieldnames(Perovskite)
@enter (1, 3) .+ (0, 7)