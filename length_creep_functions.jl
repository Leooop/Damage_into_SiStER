##########################
#### TYPES DEFINITION ####
##########################

# Parameters type
struct RockParams
    E::Float64# Young Modulus (GPa)
    μ::Float64 # Friction coef
    β::Float64 # Correction factor
    K₁c::Float64 # Critical stress intensity factor (Pa.m^(1/2))
    a::Float64 # Initial flaw size (m)
    ψ::Float64 # crack angle to the principal stress (radians)
    ρ₀::Float64 # Initial flaw density
    n::Float64 # Stress corrosion index
    n_f::Float64 # Stress corrosion index from fit
    l̇₀::Float64 # Ref. crack growth rate (m/s)
    l̇₀_f::Float64 # Ref. crack growth rate from fit (m/s)
    H::Float64# Activation enthalpy (J/mol)
    A::Float64# Preexponential factor (m/s)
    A₁::Float64 # Ashby and Sammis constant 1
    A₃::Float64 # Ashby and Sammis constant 1
end

compute_A1(β, μ) = π*sqrt(β/3)*(sqrt(1 + μ^2) + μ)
compute_A3(A1, μ) = A1*((sqrt(1 + μ^2) + μ)/(sqrt(1 + μ^2) - μ))

# Convenience functions for type instantiation :
function RockParams(E,μ,β,K₁c,a,ψ,ρ₀,n,n_f,l̇₀,l̇₀_f,H,A)
    A₁ = compute_A1(β, μ)
    A₃ = compute_A3(A₁, μ)
    RockParams(E,μ,β,K₁c,a,ψ,ρ₀,n,n_f,l̇₀,l̇₀_f,H,A,A₁,A₃)
end

RockParams(;E = 0.0,
            μ = 0.0,
            β = 0.0,
            K₁c = 0.0,
            a = 0.0,
            ψ = 0.0,
            ρ₀ = 0.0,
            n = 0.0,
            n_f = 0.0,
            l̇₀ = 0.0,
            l̇₀_f = 0.0,
            H = 0.0,
            A = 0.0
            ) = RockParams(E,μ,β,K₁c,a,ψ,ρ₀,n,n_f,l̇₀,l̇₀_f,H,A)

# Solution type :
mutable struct Solution
    t_vec::Vector{Float64}
    KI::Vector{Float64}
    l::Vector{Float64}
    ϵ₁::Vector{Float64}
    ϵ̇₁::Vector{Float64}
end

# Convenience functions :
Solution(N::Int) = Solution(Vector{Float64}(undef,N),
                            Vector{Float64}(undef,N),
                            Vector{Float64}(undef,N),
                            Vector{Float64}(undef,N),
                            Vector{Float64}(undef,N)
                            )

Solution() = Solution(Vector{Float64}(undef,0),
                            Vector{Float64}(undef,0),
                            Vector{Float64}(undef,0),
                            Vector{Float64}(undef,0),
                            Vector{Float64}(undef,0)
                            )

##############################
#### FUNCTIONS DEFINITION ####
##############################

function compute_σ₁c(p::RockParams,σ₃)
    μ, a, K₁c = p.μ, p.a, p.K₁c
    term1 = ((sqrt(1 + μ^2) + μ)/(sqrt(1 + μ^2) - μ))*σ₃
    term2 = (sqrt(3)/(sqrt(1 + μ^2) - μ))*(K₁c/sqrt(π*a))
    return term1 + term2
end

compute_Fw(p::RockParams,σ₁,σ₃) = (p.A₁*σ₁ - p.A₃*σ₃)*p.a^2

function compute_KI(p::RockParams,σ₁,σ₃,l)
    c1 = float(π)^(-2)*((l/p.a)+p.β)^(-3/2)
    c2 = (2*(π*cos(p.ψ))^(-2) * (l/p.a)^(1/2) ) / (p.ρ₀^(-2/3) - (1 + (l/(cos(p.ψ)*p.a)))^2)
    c3 = (2/π)*(l/p.a)^(1/2)
    return ((p.A₁*σ₁ - p.A₃*σ₃)*(c1 + c2) - σ₃*c3) * sqrt(π*p.a)
end

compute_c1(p::RockParams,l) = float(π)^(-2)*((l/p.a)+p.β)^(-3/2)
compute_c2(p::RockParams,l) = ( 2*(π*cos(p.ψ))^(-2) * (l/p.a)^(1/2) ) / (p.ρ₀^(-2/3) - (1 + (l/(cos(p.ψ)*p.a)))^2)
compute_c3(p::RockParams,l) = (2/π)*(l/p.a)^(1/2)

function compute_KI(p::RockParams,σ₁,σ₃,c1,c2,c3,l)
    return ((p.A₁*σ₁ - p.A₃*σ₃)*(c1 + c2) - σ₃*c3) * sqrt(π*p.a)
end

function compute_σ₁(p::RockParams,l,σ₃)
    c1 = float(π)^(-2)*((l/p.a)+p.β)^(-3/2)
    c2 = ( 2*(π*cos(p.ψ))^(-2) * (l/p.a)^(1/2) ) / (p.ρ₀^(-2/3) - (1 + (l/(cos(p.ψ)*p.a)))^2)
    c3 = (2/π)*(l/p.a)^(1/2)

    num = (p.K₁c/sqrt(π*p.a)) + σ₃*(c3 + p.A₃*(c1 + c2))
    den = p.A₁*(c1 + c2)
    return num/den
end

function compute_l0(p::RockParams,σ₁,σ₃ ; precision = 1e-6, length_max = 1e-2)
    # Make sure keyword arguments are acceptables :
    @assert (precision <= length_max) "precision must be lower than length_max"

    # synthetic crack length vector :
    l_vec = 0:precision:length_max

    # principal stress associated to it :
    sigma1_vec = compute_σ₁.(Ref(p),l_vec,σ₃)

    # make sure we are in a physically possible stress domain
    @assert (σ₁ <= maximum(sigma1_vec)) "input σ₁ is too large for the computation of an initial crack length. Maximum σ₁ should be $(maximum(sigma1_vec))"

    # index of maximum stress :
    val, ind1 = findmax(sigma1_vec)

    # index of minimum difference between σ₁ and sigma1_vec before max stress :
    val, ind2 = findmin(abs.(sigma1_vec[1:ind1] .- σ₁))

    return l_vec[ind2]
end


function compute_subcrit_growth_rate(p::RockParams, KI ; parameters = "Atkinson")
    if parameters == "Atkinson"
        return p.l̇₀*(KI/p.K₁c)^(p.n)
    elseif parameters =="data fit"
        return p.l̇₀_f*(KI/p.K₁c)^(p.n_f)
    end
end

function compute_c1c2c3(p::RockParams,l)
    c1 = float(π)^(-2)*((l/p.a)+p.β)^(-3/2)
    c2 = ( 2*(π*cos(p.ψ))^(-2) * (l/p.a)^(1/2) ) / (p.ρ₀^(-2/3) - (1 + (l/(cos(p.ψ)*p.a)))^2)
    c3 = (2/π)*(l/p.a)^(1/2)
    return c1, c2, c3
end

function iterate_adaptative_dt(p::RockParams, sol::Solution, σ₁, σ₃, Δt ; parameters = "Atkinson", e₀ = 1e-9)
    ####### computation of crack length increment over Δt using time step Δt or Δt/2 :

    # compute constants and KI at l:
    c1, c2, c3 = compute_c1c2c3(p::RockParams,sol.l[end])
    KI = compute_KI(p,σ₁,σ₃,c1,c2,c3,sol.l[end])

    KI <= 0 && return nothing # "nothing" is used to stop iterating when found in the loop in solve_strain_adaptative function

    # compute time derivative of crack length at l :
    dldt1 = compute_subcrit_growth_rate(p, KI, parameters = parameters)

    # crack growth increment over Δt :
    dl1 = dldt1*Δt

    # crack groth increment Δt/2 :
    dl_int = dldt1*(Δt/2)

    # constant and KI at l + dl_int :
    c1_int, c2_int, c3_int = compute_c1c2c3(p::RockParams,sol.l[end]+dl_int)
    KI_int = compute_KI(p,σ₁,σ₃,c1_int,c2_int,c3_int,sol.l[end]+dl_int)

    KI_int <= 0 && return nothing # "nothing" is used to stop iterating when found in the loop in solve_strain_adaptative function

    # compute time derivative of crack length at l+dl_int:
    dldt_int = compute_subcrit_growth_rate(p, KI_int, parameters = parameters)

    # crack groth increment over two times Δt/2 :
    dl2 = dl_int + dldt_int*(Δt/2)

    ######## error :
    e = abs(dl1 - dl2)
    #print(e < e₀,"\n") #Debugger print
    if e < e₀
        # compute parts of dϵ₁dt :
        term1 = (3*p.ρ₀)/(cos(p.ψ)*p.a)^3
        term2 = (sol.l[end] + cos(p.ψ)*p.a) * ((KI*sqrt(π*p.a))/p.E)
        term3 = (p.A₁*(c1 + c2))*dldt1

        l =  sol.l[end] + dl1
        dϵ₁dt = term1 * term2 * term3
        ϵ₁ = sol.ϵ₁[end] + dϵ₁dt*Δt
        t = sol.t_vec[end] + Δt
        Δt_next = min(Δt*abs(e₀/e),Δt*2)
        return t, KI, l, ϵ₁, dϵ₁dt, Δt_next
    else
        iterate_adaptative_dt(p::RockParams, sol::Solution, σ₁, σ₃, float(Δt*abs(e₀/e)^2), parameters = parameters, e₀ = e₀)
    end

end

function solve_strain_adaptative(p::RockParams, σ₁, σ₃ ; simulation_time = 5000, Δt = 1, parameters = "Atkinson", precision_l0 = 1e-6, length_max = 1e-2, abs_tol = 1e9)

    # Are we in the appropriate stress domain (ie. principal stress higher than critical stress needed to open wing cracks)?
    σ₁c = compute_σ₁c(p::RockParams,σ₃)
    σ₁ <= σ₁c && throw(DomainError("σ₁ lower than critical stress needed to open wing cracks : σ₁c = $σ₁c"))

    # Initialization of solution structure :
    sol = Solution() # Empty vector fields

    # Initial conditions :
    push!(sol.t_vec, 0.0)
    push!(sol.l, compute_l0(p,σ₁,σ₃, precision = precision_l0, length_max = length_max))
    push!(sol.KI, compute_KI(p,σ₁,σ₃,sol.l[1]))
    push!(sol.ϵ₁, 0.0)
    push!(sol.ϵ̇₁, 0.0)

    # Iteration :
    while sol.t_vec[end] <= simulation_time
        items = iterate_adaptative_dt(p::RockParams, sol::Solution, σ₁, σ₃, Δt ; parameters = parameters, e₀ = abs_tol)

        # get out of the loop if KI
        if items == nothing
            print("Simulation has stopped before reaching maximum time. Failure happened just after the last recorded iteration \n")
            break
        end
        push!(sol.t_vec, items[1])
        push!(sol.KI, items[2])
        push!(sol.l, items[3])
        push!(sol.ϵ₁, items[4])
        push!(sol.ϵ̇₁, items[5])
        Δt = items[6]
    end
    return sol
end

#=

compute_Nv(p::RockParams) = (3p.ρ₀)/(4*π*(cos(p.ψ)*p.a)^3)

function compute_σ₃i(p::RockParams,σ₁,σ₃,l)
    Fw = compute_Fw(p,σ₁,σ₃)
    S = π^(1/3)*(3/(4Nv))^(2/3)
    return Fw / (S - π*(l+cos(p.ψ)*p.a)^2)
end


function compute_max_l_increment(p::RockParams,σ₁,σ₃ ; precision = 1e-6, length_max = 1e-2, frac = 1e-3)
    # Make sure keyword arguments are acceptables :
    @assert (precision <= length_max) "precision must be lower than length_max"

    # synthetic crack length vector :
    l_vec = 0:precision:length_max

    # principal stress associated to it :
    sigma1_vec = compute_σ₁.(Ref(p),l_vec,σ₃)

    # index of maximum stress :
    val, ind = findmax(sigma1_vec)

    # critical length for macroscopic failure :
    l_max = l_vec[ind]

    # maximum length increment at each time step :
    # index of maximum stress :
    return frac * l_max
end

function compute_dϵ₁dt(p::RockParams, σ₁, σ₃, l ; parameters = "Atkinson")
    # compute constants :
    c1, c2, c3 = compute_c1c2c3(p::RockParams,l)
    # compute KI :
    KI = compute_KI(p,σ₁,σ₃,c1,c2,c3,l)

    KI <= 0 && return nothing
    # KI should be positive :
    ###KI < 0 && throw(DomainError("KI negative that should be positive. Increasing ###differential stress can help, as well as keeping σ₃ low"))

    # compute time derivative of crack length :
    dldt = compute_subcrit_growth_rate(p, KI, parameters = parameters)

    # compute parts of the solution :
    term1 = (3*p.ρ₀)/(cos(p.ψ)*p.a)^3
    term2 = (l + cos(p.ψ)*p.a) * ((KI*sqrt(π*p.a))/p.E)
    term3 = (p.A₁*(c1 + c2))*dldt
    dϵ₁dt = term1 * term2 * term3

    return KI, dldt, dϵ₁dt
end


function solve_strain(p::RockParams, σ₁, σ₃ ; N = 400000, Δt = 1, parameters = "Atkinson", precision_l0 = 1e-6, length_max = 1e-2)

    # Are we in the appropriate stress domain (ie. principal stress higher than critical stress needed to open wing cracks)?
    σ₁c = compute_σ₁c(p::RockParams,σ₃)
    σ₁ <= σ₁c && throw(DomainError("σ₁ lower than critical stress needed to open wing cracks : σ₁c = $σ₁c"))

    # Initialization of solution structure :
    sol = Solution(N)

    # Initial conditions :
    sol.t_vec[1] = 0.0
    sol.l[1] = compute_l0(p,σ₁,σ₃, precision = precision_l0, length_max = length_max)

    sol.KI[1] = compute_KI(p,σ₁,σ₃,sol.l[1])
    sol.ϵ₁[1] = 0.0
    sol.ϵ̇₁[1] = 0.0

    # Iteration :
    for i = 2:N
        KI, dldt, dϵ₁dt = compute_dϵ₁dt(p::RockParams, σ₁, σ₃, sol.l[i-1], parameters = parameters)

        sol.l[i] = sol.l[i-1] + dldt*Δt
        sol.KI[i] = KI
        sol.ϵ₁[i] = sol.ϵ₁[i-1] + dϵ₁dt*Δt
        sol.ϵ̇₁[i] = dϵ₁dt

    end
    return sol
end

function solve_strain2(p::RockParams, σ₁, σ₃ ; simulation_time = 10000, Δt = 1, parameters = "Atkinson", precision_l0 = 1e-6, length_max = 1e-2)

    # Are we in the appropriate stress domain (ie. principal stress higher than critical stress needed to open wing cracks)?
    σ₁c = compute_σ₁c(p::RockParams,σ₃)
    σ₁ <= σ₁c && throw(DomainError("σ₁ lower than critical stress needed to open wing cracks : σ₁c = $σ₁c"))

    # Initialization of solution structure :
    sol = Solution()

    # Initial conditions :
    push!(sol.t_vec, 0.0)
    push!(sol.l, compute_l0(p,σ₁,σ₃, precision = precision_l0, length_max = length_max))
    push!(sol.KI, compute_KI(p,σ₁,σ₃,sol.l[1]))
    push!(sol.ϵ₁, 0.0)
    push!(sol.ϵ̇₁, 0.0)

    # Iteration :
    while  sol.t_vec[end] <= simulation_time
        items = compute_dϵ₁dt(p::RockParams, σ₁, σ₃, sol.l[end], parameters = parameters)

        if items == nothing
            Juno.info("Simulation has stopped before reaching maximum time. Failure happened just after the last iteration recorded")
            break
        end

        push!(sol.t_vec, sol.t_vec[end]+Δt)
        push!(sol.KI, items[1])
        push!(sol.l, sol.l[end] + items[2]*Δt)
        push!(sol.ϵ₁, sol.ϵ₁[end] + items[3]*Δt)
        push!(sol.ϵ̇₁, items[3])

    end
    return sol
end

## Extention of Base.push! to operate on Solution Type :
function push!(s::Solution,items::Tuple{Float64})
    fields = propertynames(s)
    for i in length(fields)
        setproperty!(s,fields[i],items[i])
    end
end
=#
