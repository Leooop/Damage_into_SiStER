##########################
#### TYPES DEFINITION ####
##########################
# Parameters type

struct RockParams
    E::Float64 # Young Modulus (Pa)
    ν::Float64 # Poisson ratio
    λ::Float64 # first Lamé parameter (Pa)
    G::Float64 # Shear modulus (Pa)
    μ::Float64 # Friction coef
    β::Float64 # Correction factor
    K₁c::Float64 # Critical stress intensity factor (Pa.m^(1/2))
    a::Float64 # Initial flaw size (m)
    ψ::Float64 # crack angle to the principal stress (radians)
    D₀::Float64 # Initial flaw density
    n::Float64 # Stress corrosion index
    n_f::Float64 # Stress corrosion index from fit
    l̇₀::Float64 # Ref. crack growth rate (m/s)
    l̇₀_f::Float64 # Ref. crack growth rate from fit (m/s)
    H::Float64 # Activation enthalpy (J/mol)
    A::Float64 # Preexponential factor (m/s)
    A₁::Float64 # Ashby and Sammis constant 1
    A₃::Float64 # Ashby and Sammis constant 1
end

compute_A1(β, μ) = π*sqrt(β/3)*(sqrt(1 + μ^2) + μ)
compute_A3(A1, μ) = A1*((sqrt(1 + μ^2) + μ)/(sqrt(1 + μ^2) - μ))

# Convenience functions for type instantiation :
function RockParams(E,ν,μ,β,K₁c,a,ψ,D₀,n,n_f,l̇₀,l̇₀_f,H,A)
    λ = (E*ν)/((1+ν)*(1-2ν))
    G = E/(2*(1+ν))
    A₁ = compute_A1(β, μ)
    A₃ = compute_A3(A₁, μ)
    RockParams(E,ν,λ,G,μ,β,K₁c,a,ψ,D₀,n,n_f,l̇₀,l̇₀_f,H,A,A₁,A₃)
end

# Allow keyword arguments instantiation :
RockParams(;E = 0.0,
            ν = 0.0,
            μ = 0.0,
            β = 0.0,
            K₁c = 0.0,
            a = 0.0,
            ψ = 0.0,
            D₀ = 0.0,
            n = 0.0,
            n_f = 0.0,
            l̇₀ = 0.0,
            l̇₀_f = 0.0,
            H = 0.0,
            A = 0.0
            ) = RockParams(E,ν,μ,β,K₁c,a,ψ,D₀,n,n_f,l̇₀,l̇₀_f,H,A)

##############################
#### FUNCTIONS DEFINITION ####
##############################

# Eq 3 Bhat2011 :
function compute_σ₁c(p::RockParams,σ₃)
    μ, a, K₁c = p.μ, p.a, p.K₁c
    term1 = ((sqrt(1 + μ^2) + μ)/(sqrt(1 + μ^2) - μ))*σ₃
    term2 = (sqrt(3)/(sqrt(1 + μ^2) - μ))*(K₁c/sqrt(π*a))
    σ₁c = term1 - term2
    return σ₁c
end


# eq 11 Bhat2011 :
compute_Fw(p::RockParams,σ₁,σ₃) = -(p.A₁*σ₁ - p.A₃*σ₃)*p.a^2

# eq 16 Bhat2012 & 2016 & notes (because in Bhat2011 c2 isn't the same form as in Harsha's notes) :
compute_c1(p::RockParams,D) = sqrt(1-cos(p.ψ)^2)/(π*cos(p.ψ)^(3/2)*(float((D/p.D₀))^(1/3) - 1 + p.β/cos(p.ψ))^(3/2))

# Perol&Bhat2016 : 1/α  or  Harsha's notes : 1/α^2 ???
compute_c2(p::RockParams,D) = (sqrt(1 - cos(p.ψ)^2)/cos(p.ψ)^2) * (p.D₀^(2/3)/(1 - D^(2/3)))

compute_c3(p::RockParams,D) = (2/π)*sqrt(cos(p.ψ))*((D/p.D₀)^(1/3) - 1)^(1/2)

# eq 15 Bhat2012 (A1 : *c2*c3), Perol&Bhat2016 (A1 : ...*c2)*c3):
# Perol&Bhat2016 is the corrected version, and the one implemented
compute_A(p::RockParams,c1,c2,c3) = p.μ*c1 + (1 + p.μ*c2)*c3
compute_B(c1,c2,c3) = c1 + c2*c3

# eq 11 in Harsha's notes :
compute_A1(p::RockParams,A) = A * sqrt((π*p.D₀*(1 - p.ν))/cos(p.ψ)^3)
compute_B1(p::RockParams,B) = B * sqrt((π*p.D₀*(1 - p.ν))/cos(p.ψ)^3)

# Harsha's notes :
function compute_Γ(p::RockParams,A1,B1)
    term1 = (3*(1-2p.ν))/(2*(1+p.ν))
    term2 = (3*(1-2p.ν)*B1^2)/(4*(1+p.ν))
    term3 = A1^2/2
    return term1 + term2 + term3
end

function compute_c1c2c3(p::RockParams,D)
    c1 = compute_c1(p,D)
    c2 = compute_c2(p,D)
    c3 = compute_c3(p,D)
    return c1, c2, c3
end

function compute_AB(p::RockParams,c1,c2,c3)
    A = compute_A(p,c1,c2,c3)
    B = compute_B(c1,c2,c3)
    return A, B
end
function compute_AB(p::RockParams, D)
    c1, c2, c3 = compute_c1c2c3(p,D)
    A = p.μ*c1 + (1 + p.μ*c2)*c3
    B = c1 + c2*c3
    return A, B
end

function compute_A1B1(p::RockParams,A,B)
    A1 = compute_A1(p,A)
    B1 = compute_B1(p,B)
    return A1, B1
end


function compute_KI(p::RockParams,σ,τ,A,B)
    return (A*σ + B*τ) * sqrt(π*p.a)
end

function compute_KI(p::RockParams,σ,τ,D)
    c1, c2, c3 = compute_c1c2c3(p,D)
    A, B = compute_AB(p,c1,c2,c3)
    return (A*σ + B*τ) * sqrt(π*p.a)
end

#####################################################
# Computing tau max at sigma fixed isn't realistic, but this can help finding D0 a parameter search (see determine_D₀ function)
####################################################

# equating KI to KIC :
function compute_τ(p::RockParams,D,σ)
    A, B = compute_AB(p,compute_c1c2c3(p,D)...)
    τ = (1/B) * (p.K₁c/sqrt(π*p.a) - A*σ)
    return τ
end

function compute_max_dynamic_τ(rp::RockParams,D,σ)
    τ_vec = compute_τ(rp::RockParams,D,σ)
    return maximum(τ_vec)
end

function determine_D₀(rp::RockParams,D ; D₀_min =1e-4 , D₀_max = 0.2 ,precision = 1e-4)
    D₀_vec = collect(range(D₀_min,D₀_max, step = precision))
    τ_max_vec = Vector{Float64}(undef,length(D₀_vec))
    target_τ = 3.020e8
    σ = -2.043e8
    for i in eachindex(D₀_vec)
        # instantiation of the Rock Parameter object with proper D₀ :
        p = RockParams(E = rp.E,
                          ν = rp.ν,
                          μ = rp.μ,
                          β = rp.β,
                          K₁c = rp.K₁c,
                          a = rp.a,
                          ψ = rp.ψ,     #### ou π/4 (brantut)
                          D₀ = D₀_vec[i], # 2.8e-3;
                          n = rp.n,
                          n_f = rp.n_f,
                          l̇₀ = rp.l̇₀,
                          l̇₀_f = rp.l̇₀_f,
                          H = rp.H,
                          A = rp.A
                          )
        τ_max_vec[i] = maximum(compute_τ.(Ref(p),D[D .> D₀_vec[i]],σ))
    end
    val, ind = findmin(abs.(τ_max_vec .- target_τ))
    return D₀_vec[ind], τ_max_vec
end

# D_vec = collect(rp.D₀:1e-5:0.999)
# D0, tau_vec = determine_D₀(rp, D_vec, D₀_min =0.0001 , D₀_max = 0.9 ,precision = 1e-3)
# plot(range(0.0001,0.2,step = 1e-3),tau_vec)


function compute_dDdl(p::RockParams,D)
    (3*D^(2/3)*p.D₀^(1/3))/(cos(p.ψ)*p.a)
end

# dDdt :
function compute_subcrit_damage_rate(p::RockParams, KI, D ; parameters = "Atkinson")
    # crack length derivative of damage :
    dDdl = compute_dDdl(p,D)

    if parameters == "Atkinson"
        dldt = p.l̇₀*(KI/p.K₁c)^(p.n)
        KI > 0 && return dDdl * dldt
        KI <= 0 && return 0.0

    elseif parameters == "data fit"
        dldt = p.l̇₀_f*(KI/p.K₁c)^(p.n_f)
        KI > 0 && return dDdl * dldt
        KI <= 0 && return 0.0

    end
end

#### from harsha's notes :
function compute_dc1dD(p::RockParams,D)
    α = cos(p.ψ)
    term1 = (-sqrt(1-α^2)/(2*π*α^(3/2)*D^(2/3)*p.D₀^(1/3)))
    term2 = ((D/p.D₀)^(1/3) - 1 + p.β/α)^(-5/2)
    return term1 * term2
end

function compute_dc2dD(p::RockParams,D)
    α = cos(p.ψ)
    term1 = (2*sqrt(1-α^2)*p.D₀^(2/3))/(3*α^2*D^(1/3))
    term2 = (1 - D^(2/3))^(-2)
    return term1 * term2
end

function compute_dc3dD(p::RockParams,D)
    α = cos(p.ψ)
    term1 = (sqrt(α))/(3*π*D^(2/3)*p.D₀^(1/3))
    term2 = ((D/p.D₀)^(1/3) - 1)^(-1/2)
    return term1 * term2
end

function compute_dB1dD(p::RockParams,c2,c3,dc1dD,dc2dD,dc3dD)
    return (dc1dD + dc2dD*c3 + c2*dc3dD) * sqrt((π*p.D₀*(1 - p.ν))/cos(p.ψ)^3)
end

# time derivatives :
compute_dc1dt(dc1dD,dDdt) = dc1dD*dDdt
compute_dc2dt(dc2dD,dDdt) = dc2dD*dDdt
compute_dc3dt(dc3dD,dDdt) = dc3dD*dDdt
compute_dB1dt(dB1dD,dDdt) = dB1dD*dDdt





# function compute_a1(B1,Γ)
#     return (1/Γ)*(1 + B1^2/2)
# end
#
# function compute_b1(A1,B1,Γ)
#     return -(1/Γ)*((A1*B1)/2)
# end
#
# function compute_b2(p::RockParams,A1,Γ)
#     return (1/Γ)*(A1^2/2 + (3*(1-2p.ν))/(2*(1+p.ν)))
# end
#
# function compute_dΓdt(p::RockParams,A1,B1,dA1dt,dB1dt)
#     return ((3*(1-2p.ν))/(2*(1+p.ν)))*B1*dB1dt + A1*dA1dt
# end
#
# function compute_da1dt(B1,Γ,dB1dt,dΓdt)
#     return -(dΓdt/Γ^2)*(1 + B1^2/2) + (B1*dB1dt)/Γ
# end
#
# function compute_db1dt(A1,B1,Γ,dA1dt,dB1dt,dΓdt)
#     return (dΓdt/Γ^2)*((A1*B1)/2) - (1/2Γ)*(dA1dt*B1 + A1*dB1dt)
# end
#
# function compute_db2dt(p::RockParams,A1,Γ,dA1dt,dΓdt)
#     return -(dΓdt/Γ^2)*(A1^2/2 + (3*(1-2p.ν))/(2*(1+p.ν))) + (2/Γ)*A1*dA1dt
# end
#
# function compute_dσdt(p::RockParams,a1,b1,da1dt,db1dt,ϵ,γ,dϵdt,dγdt)
#     return p.μ * (da1dt*ϵ + a1*dϵdt + db1dt*γ + b1*dγdt)
# end
#
# function compute_dτdt(p::RockParams,b1,b2,db1dt,db2dt,ϵ,γ,dϵdt,dγdt)
#     return p.μ * (db1dt*ϵ + b1*dϵdt + db2dt*γ + b2*dγdt)
# end

function compute_shear_viscosity(p::RockParams,B1,dB1dt)
    η_damage = p.G/(2*B1*dB1dt)
    return η_damage
end

function compute_shear_viscosity(p::RockParams, σ, τ, D ; parameters = "Atkinson")

    c1, c2, c3 = compute_c1c2c3(p,D)
    A, B = compute_AB(p,c1,c2,c3)
    # Test stress domain to be in sub failure conditions :
    #τ_status = test_τ(p::RockParams,D,σ,τ,precision = 1e-5) # À  optimiser : recalcul de A et B
    #τ_status == false && return 1e-40

    ~, B1 = compute_A1B1(p,A,B)

    KI = compute_KI(p,σ,τ,A,B)
    # viscosity is zero if KI <=0 because cracks can't grow :
    #KI <= 0  && return 1e30
    KI < 1e-7  && (KI = 1e-7)

    dDdt = compute_subcrit_damage_rate(p, KI, D, parameters = parameters)

    dc1dD = compute_dc1dD(p,D)
    dc2dD = compute_dc2dD(p,D)
    dc3dD = compute_dc3dD(p,D)
    dB1dD = compute_dB1dD(p,c2,c3,dc1dD,dc2dD,dc3dD)

    dB1dt = compute_dB1dt(dB1dD,dDdt)

    η_damage = p.G/(2*B1*dB1dt)
    return η_damage
end

# function compute_shear_viscosity(p::RockParams, σ, τ::Vector, D ; parameters = "Atkinson")
#
#     c1, c2, c3 = compute_c1c2c3(p,D)
#     A, B = compute_AB(p,c1,c2,c3)
#     # Test stress domain to be in sub failure conditions :
#     #τ_status = test_τ(p::RockParams,D,σ,τ,precision = 1e-5) # À  optimiser : recalcul de A et B
#     #τ_status == false && return 1e-30
#
#     ~, B1 = compute_A1B1(p,A,B)
#
#     KI = compute_KI(p,σ,τ,A,B)
#     # viscosity is zero if KI <=0 because cracks can't grow :
#     #KI <= 0  && return 1e30
#     KI < 1e-7  && KI = 1e-7
#
#     dDdt = compute_subcrit_damage_rate(p, KI, D, parameters = parameters)
#
#     dc1dD = compute_dc1dD(p,D)
#     dc2dD = compute_dc2dD(p,D)
#     dc3dD = compute_dc3dD(p,D)
#     dB1dD = compute_dB1dD(p,c2,c3,dc1dD,dc2dD,dc3dD)
#
#     dB1dt = compute_dB1dt(dB1dD,dDdt)
#
#     η_damage = p.G/(2*B1*dB1dt)
#     return η_damage
# end

function compute_γ̇_viscosity(rp::RockParams,σ::Float64,τ::Float64,D::Float64 ; parameters = "Atkinson")
    η = compute_shear_viscosity(rp, σ, τ, D, parameters = parameters)
    γ̇ = τ./(2η)
    return γ̇
end


# function compute_γ̇_viscosity(rp::RockParams,σ::Float64,τ::Vector,D::Float64 ; parameters = "Atkinson")
#     η = compute_shear_viscosity.(Ref(rp), σ::Float64, τ, D, parameters = parameters)
#     γ̇ = τ./(2η)
#     return γ̇
# end

function compute_γ̇_viscosity(rp::RockParams,σ₃::Vector,σ₁::Vector,D; parameters = "Atkinson")
    γ̇_mat = Matrix{Float64}(undef,length(σ₁),length(σ₃))
    σ_mat = Matrix{Float64}(undef,length(σ₁),length(σ₃))
    τ_mat = Matrix{Float64}(undef,length(σ₁),length(σ₃))
    for j in eachindex(σ₃), i in eachindex(σ₁)
        σ_mat[i,j], τ_mat[i,j] = stress2invariants(σ₁[i],σ₃[j])
        γ̇_mat[i,j] = compute_γ̇_viscosity(rp, σ_mat[i,j], τ_mat[i,j], D, parameters = parameters)
    end
    return σ_mat, τ_mat, γ̇_mat
end

function compute_γ̇_viscosity(rp::RockParams,σ₃::Float64,σ₁::Vector,D::Vector; parameters = "Atkinson")
    γ̇_mat = Matrix{Float64}(undef,length(σ₁),length(D))
    σ_mat = Matrix{Float64}(undef,length(σ₁),length(D))
    τ_mat = Matrix{Float64}(undef,length(σ₁),length(D))
    for j in eachindex(D), i in eachindex(σ₁)
        σ_mat[i,j], τ_mat[i,j] = stress2invariants(σ₁[i],σ₃)
        γ̇_mat[i,j] = compute_γ̇_viscosity(rp, σ_mat[i,j], τ_mat[i,j], D[j], parameters = parameters)
    end
    return σ_mat, τ_mat, γ̇_mat
end



function plot_γ̇_viscosity(σ₃::Vector,τ::Matrix,γ̇::Matrix,D::Float64 ; log_τ = true)
    # Careful \sigma3 is constant over each row computed, but σ isn't.
    # plots strain rate as a function of differential stress, not strain invariant

    # Choice of expression of tau :
    τ = τ.*sqrt(3) # converts to differential stress
    log_τ && (τ = log10.(τ))
    # Color palette :
    colors = [:blue,:red,:green,:magenta,:orange,:yellow,:black]

    # Plots :
    p = plot(τ[:,1],log10.(γ̇[:,1]),label = L"σ₃="*"$(σ₃[1])", c=colors[1])
    for j = 2:size(γ̇,2)
        plot!(τ[:,j],log10.(γ̇[:,j]),label = L"σ₃="*"$(σ₃[j])", c=colors[j])
    end

    # Labels :
    log_τ ? xlabel!(L"log(Δσ) (Pa)") : xlabel!(L"Δσ (Pa)")
    ylabel!(L"log(γ̇)")


    # plot brittle creep data :
    Δσ_bc = [392.2,399.5,408.5,424.3,434.5,450.9,466.2,484.9,494.7].*1e6
    σ₃_bc = zeros(length(Δσ_bc)) .- 30e6
    ϵ̇_bc = [5.33e-9,6.91e-9,2.34e-8,4.27e-8,9.89e-8,1.73e-7,2.20e-7,3.81e-7,8.14e-7]

    # plot points :
    log_τ && (Δσ_bc = log10.(Δσ_bc))
    scatter!(Δσ_bc,log10.(ϵ̇_bc),color = :red, label = "brittle creep data")

    # Plot max τ :
    # y_vec = collect(range(1e-30,1e40,length = 11))
    # for j = 1:size(γ̇,2)
    #     max_τ = compute_max_dynamic_τ(rp,D,σ[j])
    #     max_τ_vec = zeros(length(y_vec)) .+ max_τ
    #     log_τ && (max_τ_vec = log10.(max_τ_vec))
    #     plot!(max_τ_vec,log10.(y_vec),label="",c=colors[j],ls = :dash)
    # end
    return p
end

function plot_γ̇_viscosity(σ₃::Float64,τ::Matrix,γ̇::Matrix,D::Vector ; log_τ = true)
    # Choice of expression of tau :
    τ = τ.*sqrt(3)
    log_τ && (τ = log10.(τ))
    # Color palette :
    colors = [:blue,:red,:green,:magenta,:orange,:yellow,:black]

    # Plots :
    gamma_dot = γ̇[:,1]
    tau = τ[:,1]
    p = plot(tau[γ̇[:,1].>=0],log10.(gamma_dot[γ̇[:,1].>=0]),label = L"D="*@sprintf("%.2f",D[1]), c=colors[1])
    for j = 2:size(γ̇,2)
        gamma_dot = γ̇[:,j]
        tau = τ[:,j]
        #plot!(τ[:,j].*sqrt(3),log10.(γ̇[:,j]),label = L"D="*"$(D[j])", c=colors[j])
        plot!(tau[γ̇[:,j].>=0],log10.(gamma_dot[γ̇[:,j].>=0]),label = L"D="*@sprintf("%.2f",D[j]), c=colors[j])
    end

    # Labels :
    log_τ ? xlabel!(L"log(Δσ) (Pa)") : xlabel!(L"Δσ (Pa)")
    ylabel!(L"log(γ̇)")

    # plot brittle creep data :
    Δσ_bc = [392.2,399.5,408.5,424.3,434.5,450.9,466.2,484.9,494.7].*1e6
    σ₃_bc = zeros(length(Δσ_bc)) .- 30e6
    ϵ̇_bc = [5.33e-9,6.91e-9,2.34e-8,4.27e-8,9.89e-8,1.73e-7,2.20e-7,3.81e-7,8.14e-7]

    # plot points :
    log_τ && (Δσ_bc = log10.(Δσ_bc))
    scatter!(Δσ_bc,log10.(ϵ̇_bc),color = :red, label = "brittle creep data")

    # Plot max τ :
    # y_vec = collect(range(1e-30,1e40,length = 11))
    # for j = 1:size(γ̇,2)
    #     max_τ = compute_max_dynamic_τ(rp,D[j],σ)
    #     max_τ_vec = zeros(length(y_vec)) .+ max_τ
    #     log_τ && (max_τ_vec = log10.(max_τ_vec))
    #     plot!(max_τ_vec,log10.(y_vec),label="",c=colors[j],ls = :dash)
    # end
    # return p
end

rp = RockParams(E = 70e9,
              ν = 0.3,
              μ = 0.64, # 0.64
              β = 0.1, # 1.0
              K₁c = 1e6,
              a = 1e-3, #5e-3
              ψ = 0.5*atan(1/0.64),     #### ou π/4, 0.5*atan(1/0.64) == deg2rad(45 - 0.5*atand(0.64))
              D₀ = 0.1621, # 2.8e-3;
              n = 12,
              n_f = 5.5,
              l̇₀ = 100e-3, #100e-3
              l̇₀_f = 7e-5*10000,
              H = 50e3,
              A = 5.71
              )

# Parameters :
σ₃ = -30e6
σ₁_vec = collect(range(σ₃,-8e8,length = 50))
σ₃_vec = [0,-5e6,-10e6,-15e6,-20e6,-25e6,-30e6]
D = 0.5
D_vec = collect(range(rp.D₀+0.001,0.9,length = 7))

# at D const :
σ_mat, τ_mat, γ̇_mat = compute_γ̇_viscosity(rp,σ₃_vec,σ₁_vec,D, parameters = "Atkinson")
p = plot_γ̇_viscosity(σ₃_vec,τ_mat,γ̇_mat,D,log_τ=false)

# at σ const :
σ_mat, τ_mat, γ̇_mat = compute_γ̇_viscosity(rp,σ₃,σ₁_vec,D_vec, parameters = "data fit")
p = plot_γ̇_viscosity(σ₃,τ_mat,γ̇_mat,D_vec,log_τ=false)



#compute_γ̇_viscosity(rp,-160e6,110e6,1e-3 ; parameters = "data fit")

#KI = compute_KI(rp,-110e6,0.64*110e6,rp.D₀)
