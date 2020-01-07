"""

    fdt(λ::Real, ν::Real, dt::Real)

Compute the variance that satisfies the fluctuation-dissipation theorem.

# External links

* [Fluctuation-dissipation theorem on Wikipedia](https://en.wikipedia.org/wiki/Fluctuation-dissipation_theorem)

# Examples
```jldoctest
julia> fdt(1.0, 1.0, 0.1)
1.378404875209022

```

"""
function fdt(λ::Real, ν::Real, dt::Real)
    return sqrt(ν * (2*λ*dt - λ^2*dt^2) / dt)
end


"""

    CAM1D(d::SGS, n::Integer; dt::Real=1/24, λ::Real=1.0)

Create a timeseries using the one-dimensional CAM noise model of Sardeshmukh and
Sura (2009) for a particular SGS distribution.

See also:

>Sardeshmukh, P. D., and P. Sura, 2009: Reconciling Non-Gaussian Climate Statistics with Linear Dynamics. *Journal of Climate*, **22**, 1193–1207, [https://doi.org/10.1175/2008JCLI2358.1](https://doi.org/10.1175/2008JCLI2358.1).
>
>Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. *Journal of Climate*, **28**, 9166–9187, [https://doi.org/10.1175/JCLI-D-15-0020.1](https://doi.org/10.1175/JCLI-D-15-0020.1).

# Examples
```jldoctest
julia> d = fit(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = CAM1D(d, 1000, seed=42)
retcode: Success
Interpolation: 1st order linear
t: 1001-element Array{Float64,1}:
  0.0
  0.041666666666666664
  0.08333333333333333
  0.125
  0.16666666666666666
  0.20833333333333331
  0.24999999999999997
  0.29166666666666663
  0.3333333333333333
  0.375
  ⋮
 41.333333333333165
 41.37499999999983
 41.416666666666494
 41.45833333333316
 41.49999999999982
 41.54166666666649
 41.58333333333315
 41.624999999999815
 41.625
u: 1001-element Array{Array{Float64,1},1}:
 [0.0]
 [-0.05538507236046931]
 [0.3759078513765659]
 [0.31482891888314446]
 [0.2586995978771358]
 [0.0063144005340343146]
 [0.0819076515551697]
 [0.07880240481554046]
 [-0.09903540075387843]
 [-0.0867782526753085]
 ⋮
 [-0.048583800782025996]
 [0.30770620643597935]
 [0.6914165604620269]
 [1.2048828846162343]
 [0.6995798768325538]
 [0.5417553247799438]
 [0.593919650292867]
 [0.6785936302181321]
 [0.6785940957553152]

```

"""
function CAM1D(d::SGS, n::Integer; dt::Real=1/24, λ::Real=1.0, seed::Real=-1)
    # TODO: What are the best values for the default arguments?

    tspan = (0.0, n*dt-dt)      # time domain
    u0 = [0.0]                  # initial value
    p = [d.E, d.b, d.g, λ]      # parameter array

    prob = SDEProblem(CAM1D_dt, CAM1D_dW, u0, tspan, p, noise_rate_prototype=zeros(1,2))  #,reltol=1e-8,abstol=1e-8)

    # alg = RKMil(interpretation=:Stratonovich)
    alg = EulerHeun()

    if seed > -1
        sol = solve(prob, alg, dt=dt, seed=seed, abstol=1e-10, reltol=1e-10)
    else
        sol = solve(prob, alg, dt=dt, abstol=1e-10, reltol=1e-10) #, seed=Random.seed!())
    end

    return sol

end


function CAM1D(u0::AbstractArray, tspan::Tuple, p::AbstractArray, dt::Real; seed::Real=-1)

    prob = SDEProblem(CAM1D_dt, CAM1D_dW, u0, tspan, p, noise_rate_prototype=zeros(1,2))  #,reltol=1e-8,abstol=1e-8)

    if seed > -1
        sol = solve(prob, EulerHeun(), dt=dt, seed=seed)
    else
        sol = solve(prob, EulerHeun(), dt=dt) #, seed=Random.seed!())
    end

    return sol

end


"""

    Hasselmann1D(n::Integer; dt::Real=1/24, λ::Real=1.0, seed::Real=-1)

Create a timeseries using the one-dimensional noise model of Hasselmann (1976).

See also:

>Hasselmann, K., 1976: Stochastic climate models Part I. Theory. *Tellus*, **28**, 473–485, [https://doi.org/10.1111/j.2153-3490.1976.tb00696.x](https://doi.org/10.1111/j.2153-3490.1976.tb00696.x).

# Examples
```jldoctest
julia> Hasselmann1D(1000, seed=42)
retcode: Success
Interpolation: 1st order linear
t: 1001-element Array{Float64,1}:
  0.0
  0.041666666666666664
  0.08333333333333333
  0.125
  0.16666666666666666
  0.20833333333333331
  0.24999999999999997
  0.29166666666666663
  0.3333333333333333
  0.375
  ⋮
 41.333333333333165
 41.37499999999983
 41.416666666666494
 41.45833333333316
 41.49999999999982
 41.54166666666649
 41.58333333333315
 41.624999999999815
 41.625
u: 1001-element Array{Array{Float64,1},1}:
 [0.0]
 [-0.13005535542166108]
 [0.04828669603525988]
 [0.658691166050971]
 [0.4766431404349038]
 [0.2527905702593073]
 [0.5063085613562422]
 [0.4424101651546515]
 [0.42457583766050366]
 [0.17019337828909467]
 ⋮
 [0.04255470162528974]
 [-0.009396010196968833]
 [-0.4700188956712661]
 [-0.6640224731670598]
 [-0.6412697125734018]
 [-0.8226707386234939]
 [-0.9240097641730844]
 [-0.3550868399241892]
 [-0.35508653697280324]

```

"""
function Hasselmann1D(n::Integer; dt::Real=1/24, λ::Real=1.0, seed::Real=-1)
    # TODO: What are the best values for the default arguments?

    ν = 1.0                     # ASSUMED: variance of the white noise is 1.0
    σ = fdt(λ, ν, dt)           # scale the variance for the Hasselman case.
    p = [λ, σ]                  # DiffEq.jl parameter array
    tspan = (0.0, n*dt-dt)      # time domain
    u0 = [0.0]                  # initial value

    prob = SDEProblem(Hasselmann1D_dt, Hasselmann1D_dW, u0, tspan, p)

    if seed > -1
        sol = solve(prob, EM(), dt=dt, seed=seed)
    else
        sol = solve(prob, EM(), dt=dt)
    end

    return sol

end


function Hasselmann1D(u0::AbstractArray, tspan::Tuple, p::AbstractArray, dt::Real; seed::Real=-1)

    prob = SDEProblem(Hasselmann1D_dt, Hasselmann1D_dW, u0, tspan, p)

    if seed > -1
        sol = solve(prob, EM(), dt=dt, seed=seed)
    else
        sol = solve(prob, EM(), dt=dt)
    end

    return sol

end


function CAM1D_dt(du::AbstractArray, u::AbstractArray, p::AbstractArray, t::Real)

    # Unpack the parameters.
    E = p[1]
    b = p[2]
    g = p[3]
    λ = p[4]

    # Unpack the variables.
    x = u[1]

    du[1] = -((1.0 + 0.5 * E^2) * x + 0.5 * E * g) * λ

    return du

end

function CAM1D_dW(du::AbstractArray, u::AbstractArray, p::AbstractArray, t::Real)

    # Unpack the parameters.
    E = p[1]
    b = p[2]
    g = p[3]
    λ = p[4]

    # Unpack the variables.
    x = u[1]

    # Calculate two temporally uncorrelated standard normal variables.
    # η1 = randn()
    # η2 = randn()

    # du[1] = (b * η1 + (E * x + g) * η2) * sqrt(λ)
    du[1,1] = b * sqrt(λ)
    du[1,2] = (E * x + g) * sqrt(λ)

    return du

end


"""

    CAM1D_dx(d::SGS, x::Real, λ::Real, η1::Float64, η2::Float64, dt::Real)

Compute the change in the stochastically perturbed damped linear Markov process
``x(t)`` in a small time interval, ``dt``, for the CAM noise model.

"""
function CAM1D_dx(d::SGS, x::Real, λ::Real, η1::Float64, η2::Float64, dt::Real)
    (E, b, g) = params(d)

    # Construct the change, dx, of the stochastically perturbed damped linear
    # Markov process x(t) in a small time interval dt.
    term1 = -((1.0 + 0.5 * E^2) * x + 0.5 * E * g) * λ * dt
    term2 = ((b * η1) + ((E * x + g) * η2)) * sqrt(λ * dt)

    return term1 + term2

end


function Hasselmann1D_dt(du::AbstractArray, u::AbstractArray, p::AbstractArray, t::Real)

    # Unpack the parameters, variables.
    λ = p[1]
    σ = p[2]
    x = u[1]

    du[1] = -(λ * x)

    return du

end


function Hasselmann1D_dW(du::AbstractArray, u::AbstractArray, p::AbstractArray, t::Real)

    # Unpack the parameters, variables.
    λ = p[1]
    σ = p[2]
    x = u[1]

    du[1] = σ

    return du

end


"""

    Hasselmann1D_dx(x::Real, λ::Real, σ::Real, η::Float64, dt::Real)

Compute the change in the stochastically perturbed damped linear Markov process
``x(t)`` in a small time interval, ``dt``, for the Hasselmann model.

"""
function Hasselmann1D_dx(x::Real, λ::Real, σ::Real, η::Float64, dt::Real)

    return (-λ * x * dt) + (σ * η * sqrt(dt))

end


function Euler(d::Normal, u0::AbstractArray, tspan::Tuple, p::AbstractArray, dt::Real)

    # Unpack the parameters, variables.
    λ = p[1]
    σ = p[2]
    n = Int(round((tspan[2] - tspan[1]) / dt))     # number of iterations

    η = randn(n)
    u = Array{Float64}(undef, n)
    t = Array{Float64}(undef, n)

    u[1] = u0[1]
    t[1] = tspan[1]

    for i=2:n
        u[i] = u[i-1] + Hasselmann1D_dx(u[i-1], λ, σ, η[i], dt)
        t[i] = t[i-1] + dt
    end

    return (u, t)

end


function Euler(d::SGS, u0::AbstractArray, tspan::Tuple, p::AbstractArray, dt::Real)

    # Unpack the parameters.
    E = p[1]
    b = p[2]
    g = p[3]
    λ = p[4]
    n = Int(round((tspan[2] - tspan[1]) / dt))       # number of iterations

    # Calculate two temporally uncorrelated standard normal variables.
    η1 = randn(n)
    η2 = randn(n)
    u = Array{Float64}(undef, n)
    t = Array{Float64}(undef, n)

    u[1] = u0[1]
    t[1] = tspan[1]

    for i=2:n
        u[i] = u[i-1] + CAM1D_dx(d, u[i-1], λ, η1[i], η2[i], dt)
        t[i] = t[i-1] + dt
    end

    return (u, t)

end


function RK4(d::Normal, u0::AbstractArray, tspan::Tuple, p::AbstractArray, dt::Real)

    # Unpack the parameters, variables.
    λ = p[1]
    σ = p[2]
    n = Int(round((tspan[2] - tspan[1]) / dt))       # number of iterations

    η = randn(n)
    u = Array{Float64}(undef, n)
    t = Array{Float64}(undef, n)

    u[1] = u0[1]
    t[1] = tspan[1]

    for i=2:n
        k11 = Hasselmann1D_dx(u[i-1], λ, σ, η[i], dt)
        k21 = Hasselmann1D_dx(u[i-1] + (k11 / 2.0), λ, σ, η[i], dt)
        k31 = Hasselmann1D_dx(u[i-1] + (k21 / 2.0), λ, σ, η[i], dt)
        k41 = Hasselmann1D_dx(u[i-1] + k31, λ, 1.0, η[i], dt)

        u[i] = u[i-1] + (k11 + 2.0*k21 + 2.0*k31 + k41) / 6.0
        t[i] = t[i-1] + dt
    end

    return (u, t)

end


function RK4(d::SGS, u0::AbstractArray, tspan::Tuple, p::AbstractArray, dt::Real)

    # Unpack the parameters.
    E = p[1]
    b = p[2]
    g = p[3]
    λ = p[4]
    n = Int(round((tspan[2] - tspan[1]) / dt))       # number of iterations

    # Calculate two temporally uncorrelated standard normal variables.
    η1 = randn(n)
    η2 = randn(n)
    u = Array{Float64}(undef, n)
    t = Array{Float64}(undef, n)

    u[1] = u0[1]
    t[1] = tspan[1]

    for i=2:n
        # Do not mulitply kij by dt?
        k11 = CAM1D_dx(d, u[i-1], λ, η1[i], η2[i], dt)
        k21 = CAM1D_dx(d, u[i-1] + (k11 / 2.0), λ, η1[i], η2[i], dt)
        k31 = CAM1D_dx(d, u[i-1] + (k21 / 2.0), λ, η1[i], η2[i], dt)
        k41 = CAM1D_dx(d, u[i-1] + k31, λ, η1[i], η2[i], dt)

        u[i] = u[i-1] + (k11 + 2.0*k21 + 2.0*k31 + k41) / 6.0
        t[i] = t[i-1] + dt
    end

    return (u, t)

end
