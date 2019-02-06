## put this file's directory as working directory
using Plots; gr()


# Rock properties
wg_p = RockParams(E = 70e9,
                  μ = 0.7,
                  β = 1.15,
                  K₁c = 1.74e6,
                  a = 1.1e-3,
                  ψ = π/4,
                  ρ₀ = 2.8e-3,
                  n = 34.0,
                  n_f = 5.5,
                  l̇₀ = 0.24,
                  l̇₀_f = 7e-5,
                  H = 50e3,
                  A = 5.71
                  )

# Chosen parameters :
σ₁, σ₃ = 380e6, 30e6
max_time = 10000000
Δt_ini = 0.000001

# Solve :

sol = solve_strain_adaptative(wg_p, σ₁, σ₃ ; simulation_time = max_time, Δt = Δt_ini, parameters = "Atkinson", precision_l0 = 1e-9, length_max = 1e-2, abs_tol = 1e-9)


# plot KI/KIC

plot(sol.t_vec./3600, sol.KI./wg_p.K₁c)
    ylims!((0.5,1.5))
    xlabel!("Time (hours)")
    ylabel!("KI/KIC")

# plot l :
plot(sol.t_vec./3600, sol.l.*1000)
    xlabel!("Time (hours)")
    ylabel!("crack length (mm)")
# plot ϵ :
plot(sol.t_vec[sol.ϵ₁.<2e-5]./3600, sol.ϵ₁[sol.ϵ₁.<2e-5]) #[sol.ϵ₁.<5e-3]
    xlabel!("Time (hours)")
    ylabel!("Axial strain")

plot
sol.KI
#compute_l0(wg_props,σ₁,σ₃, precision = 1e-9, length_max = 1e-2)

minimum(sol.KI)
