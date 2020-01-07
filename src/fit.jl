"""

    fit(d::Type{SGS}, x::AbstractArray)

Fit an SGS distribution to the time series, x.

# Examples
```jldoctest
julia> fit(SGS, air())
SGS{Float64}(E=0.36623331624359673, b=0.7133796999284913, g=-1.1648873051087967)

```

"""
function fit(d::Type{SGS}, x::AbstractArray)

    # Calculate the moments.
    variance = var(x, corrected=false)
    skew = skewness(x)
    kurt = kurtosis(x)

    # Test for the skewness-kurtosis relationship constraint.
    skewbound = 1.5 * skew^2
    @debug "K-S inequality constraint: $kurt > $skewbound"
    kurt > skewbound || @warn "K-S inequality constraint violated."

    # Fit the SGS distribution using the method of moments.
    return fit_mm(d, variance, skew, kurt)

end

"""

    fit(d::Type{SGS}, variance::Real, skew::Real, kurt::Real)

Fit an SGS distribution to the desired variance, skewness and kurtosis.

# Examples
```jldoctest
julia> fit(SGS, 1.0, 1.0, 5.0)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

```

"""
function fit(d::Type{SGS}, variance::Real, skew::Real, kurt::Real)

    # Test for the skewness-kurtosis relationship constraint.
    skewbound = 1.5 * skew^2
    kurt > skewbound || @warn "K-S inequality constraint violated: $kurt > $skewbound"

    # Fit the SGS distribution using the method of moments.
    return fit_mm(d, variance, skew, kurt)

end

"""

    fit(d::Type{SGS}, variance::Real, skew::Real, kurt::Real)

Use a method of moments to fit an SGS distribution to the desired variance,
skewness and kurtosis.

# Examples
```jldoctest
julia> fit_mm(SGS, 1.0, 1.0, 5.0)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

```

"""
function fit_mm(d::Type{SGS}, variance::Real, skew::Real, kurt::Real)

    # Calculate the SGS parameters using the calculated moments.
    E = param_E(skew, kurt)
    b = param_b(variance, skew, E)
    g = param_g(variance, skew, E)

    @debug "E = $E"
    @debug "b = $b"
    @debug "g = $g"
    @debug "variance = $variance"
    @debug "skewness = $skew"
    @debug "kurtosis = $kurt"

    return d(E, b, g)

end


"""

    param_E(skew::Real, kurt::Real)

Solve for the parameter E in Equation 8(a) of Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> param_E(1.0, 5.0)
0.6236095644623235

```

"""
function param_E(skew::Real, kurt::Real)

    E = _param_E(skew, kurt)

    # TODO: Should the value of the kurtosis changed globally or just for the
    # purpose of calculating E?

    # Must implement the epsilon constraint.
    epsilon_constraint(skew, E) || @warn "Epsilon constraint for parameter E not satisfied, incrementing kurtosis"
    # epsilon_constraint(skew, E) || @debug "Epsilon constraint for parameter E not satisfied, incrementing kurtosis"
    while epsilon_constraint(skew, E) == false
        # println("epsilon constraint unmet, incrementing kurtosis")
        kurt += 0.1
        E = _param_E(skew, kurt)
        @debug "E = $E, kurt = $kurt"
        if kurt > 15
            @warn "E = $E  skew=$skew  kurt = $kurt"
            return E   # not great, but needs to be investigated
        end
    end

    return E

end

"""

    param_E(skew::AbstractArray, kurt::AbstractArray)

Solve for the parameter E in Equation 8(a) of Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> param_E([1.0, 1.2, 1.4], [3.0, 4.0, 5.0])
3-element Array{Float64,1}:
 0.5
 0.5186577368103327
 0.5220026556319158

```

"""
function param_E(skew::AbstractArray, kurt::AbstractArray)

    E = _param_E.(skew, kurt)

    # TODO: Must implement the epsilon constraint.

    return E

end

function _param_E(skew::Real, kurt::Real)

    E2 = (2.0 / 3.0) * ((kurt - (1.5 * skew^2)) / (kurt - skew^2 + 2.0))

    if E2 ≤ 0.0
        @debug "Constraint E ≥ 0 not satisfied, setting E to NaN"
        return NaN
    else
        return sqrt(E2)
    end

end

function _param_E(skew::AbstractArray, kurt::AbstractArray)

    E = sqrt.((2.0 / 3.0) * ((kurt .- (1.5 .* skew .^2)) / (kurt .- skew .^2 + 2.0)))
    return E

    # TODO: Implement the positive definite constraint on E here.

end

"""

    param_g(variance::Real, skew::Real, E::Real)

Solve for the parameter g in Equation 8(b) of Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> E = 0.6236095644623235
0.6236095644623235

julia> param_g(1.0, 1.0, E)
0.48997894350611143

```

"""
function param_g(variance::Real, skew::Real, E::Real)

    return skew * sqrt(variance) * ((1.0 - E^2) / (2.0 * E))

end

"""

    param_g(variance::AbstractArray, skew::AbstractArray, E::AbstractArray)

Solve for the parameter g in Equation 8(b) of Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> E = 0.6236095644623235
0.6236095644623235

julia> param_g(1.0, 1.0, E)
0.48997894350611143

```

"""
function param_g(variance::AbstractArray, skew::AbstractArray, E::AbstractArray)

    return skew .* sqrt.(variance) .* ((1.0 .- E .^2) / (2.0 .* E))

end

"""

    param_b(variance::Real, skew::Real, E::Real)

Solve for the parameter b in Equation 8(c) of Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> E = 0.6236095644623235
0.6236095644623235

julia> param_b(1.0, 1.0, E)
1.1709106481844573

```

"""
function param_b(variance::Real, skew::Real, E::Real)

    b2 = 2.0 * variance * (1.0 - (E^2 / 2.0) - (skew^2 * ((1.0 - E^2)^2 / (8.0 * E^2))))
    if b2 < 0.0
        @warn "Constraint b > 0 not satisfied, setting b to NaN"
        return NaN
    else
        return sqrt(b2)
    end

end

"""

    param_b(variance::AbstractArray, skew::AbstractArray, E::AbstractArray)

Solve for the parameter b in Equation 8(c) of Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> E = 0.6236095644623235
0.6236095644623235

julia> param_b(1.0, 1.0, E)
1.1709106481844573

```

"""
function param_b(variance::AbstractArray, skew::AbstractArray, E::AbstractArray)

    b2 .= 2.0 .* variance .* (1.0 .- (E .^2 / 2.0) .- (skew .^2 .* (1.0 .- E .^2) .^2 ./ (8.0 .* E .^2)))
    return b2 .< 0.0 ? NaN : sqrt.(b2)

end
