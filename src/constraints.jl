"""

    moment_constraint(order::Integer, E::Real)

Compute the moment constraint for a given order.

# Notes
Sardeshmukh et al. (2015) states that a necessary condition for the moment
``\\langle x^{n} \\rangle`` to exist (where `n` is the moment order) is that
``E^{2} < \\frac{2}{n-1}``.

"""
function moment_constraint(order::Integer, E::Real)
    return E^2 < 2 / (order - 1)
end


"""

    small_E_constraint(x::Real, skew::Real)

Compute the small E constraint for the first-order approximation to the SGS
pdf in the limit as E approaches zero.

# Examples
```jldoctest
julia> small_E_constraint(0.0, 1.0)
true

```

# Notes
Sardeshmukh et al. (2015) states that in order for the first-order approximation
of the SGS pdf to be valid in the limit ``E \\to 0``, ``|Sx| \\le 2`` where `x` is
the standardized anomaly.

"""
function small_E_constraint(x::Real, skew::Real)
    return abs(skew * x) ≤ 2
end

"""

    epsilon_condition(skew::Real)

Solve for the epsilon condition used in the epsilon constraint.

# Examples
```jldoctest
julia> epsilon_condition(1.0)
1.118033988749895

```

"""
function epsilon_condition(skew::Real)
    return sqrt(1.0 + (skew^2 / 4.0))
end

"""

    epsilon_condition(skew::AbstractArray)

Solve for the epsilon condition used in the epsilon constraint.

# Examples
```jldoctest
julia> epsilon_condition([1, 3, 4])
3-element Array{Float64,1}:
 1.118033988749895
 1.8027756377319946
 2.23606797749979

```

"""
function epsilon_condition(skew::AbstractArray)
    return sqrt.(1.0 .+ (skew .^2 ./ 4.0))
end


"""

    epsilon_constraint(skew::Real, kurt::Real)

Test the epsilon constraint used in the method of moments calculation for
parameter E in Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> epsilon_constraint(1, 0.5)
true

```

"""
function epsilon_constraint(skew::Real, E::Real)
    return E^2 > (epsilon_condition(skew) - 1.0) / epsilon_condition(skew)
end

"""

    epsilon_constraint(skew::Real, kurt::Real)

Test the epsilon constraint used in the method of moments calculation for
parameter E in Sardeshmukh et al. (2015).

# Examples
```jldoctest
julia> epsilon_constraint([1, 2, 3], [0.5, 1, 1.5])
3-element BitArray{1}:
 1
 1
 1

```

"""
function epsilon_constraint(skew::AbstractArray, E::AbstractArray)
    return E .^2 .> (epsilon_condition(skew) .- 1.0) ./ epsilon_condition(skew)
end


"""

    ks_constraint(skew::Real, kurt::Real)

Test the skewness-kurtosis inequality.

# Examples
```jldoctest
julia> ks_constraint(1, 1)
false

```

"""
function ks_constraint(skew::Real, kurt::Real)
    return kurt > 1.5 * (skew^2)
end


"""

    ks_constraint(skew::Real, kurt::Real)

Test the skewness-kurtosis inequality.

# Examples
```jldoctest
julia> ks_constraint([0.5, 1, 1.5], [3, 3, 3])
3-element BitArray{1}:
 1
 1
 0

```

"""
function ks_constraint(skew::AbstractArray, kurt::AbstractArray)
    return kurt .> 1.5 .* (skew .^2)
end
