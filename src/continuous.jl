"""

    SGS(E, b, g)

The stochastically generated skewed (SGS) distribution with parameters E, b and g, derived by Sardeshmukh and Sura (2009) and following the form of Sardeshmukh et al. (2015). The SGS probability density function (pdf) is

```math
p(x) = \\frac{1}{\\mathcal{N}} \\left[ \\left( Ex + g \\right)^{2} + b^{2}  \\right]^{-\\left(1 + \\left(1/E^{2}\\right)\\right)} \\mathrm{exp} \\! \\left[ \\frac{2g}{E^{2}b} \\mathrm{arctan} \\! \\left(\\frac{Ex + g}{b} \\right) \\right]
```

where the normalization constant, ``\\mathcal{N}``, is given as

```math
\\mathcal{N} = \\frac{2 \\pi \\nu^{1/2} \\left( 2b \\right)^{-\\left(2 \\nu + 1 \\right)} \\Gamma \\! \\left(2\\nu + 1\\right)}{\\Gamma \\! \\left( \\nu + 1 - iq/2 \\right) \\Gamma \\! \\left( \\nu + 1 + iq/2 \\right)}
```

See also:

>Sardeshmukh, P. D., and P. Sura, 2009: Reconciling Non-Gaussian Climate Statistics with Linear Dynamics. *Journal of Climate*, **22**, 1193–1207, [https://doi.org/10.1175/2008JCLI2358.1](https://doi.org/10.1175/2008JCLI2358.1).
>
>Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. *Journal of Climate*, **28**, 9166–9187, [https://doi.org/10.1175/JCLI-D-15-0020.1](https://doi.org/10.1175/JCLI-D-15-0020.1).

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

```

"""
struct SGS{T<:Real} <: ContinuousUnivariateDistribution
    E::T
    b::T
    g::T

    function SGS{T}(E::T, b::T, g::T) where {T}
        # This constraint check should be required if estimating the SGS
        # distribution using the method of moments.
        # @check_args(SGS, E ≥ zero(E) && b > zero(b))
        new(E, b, g)
    end
end

# Outer constructors.
SGS(E::T, b::T, g::T) where {T<:Real} = SGS{T}(E, b, g)
# NOTE: dots after promote()... expand the tuple.
SGS(E::Real, b::Real, g::Real) = SGS(promote(E, b, g)...)
SGS(E::Integer, b::Integer, g::Integer) = SGS(Float64(E), Float64(b), Float64(g))
SGS() = SGS(1.0, 1.0, 1.0)

# Type conversions.
# TODO: More type checking here.

# Still can't seem to find much on the @distr_support macro.
@distr_support SGS -Inf Inf

# Parameter retrieval
# See https://juliastats.github.io/Distributions.jl/stable/univariate/#Parameter-Retrieval-1

"""

    params(d::SGS)

Return a tuple of SGS parameters `E`, `b` and `g`. Let `d` be a distribution of type `D`, then `D(params(d)...)` will construct exactly the same distribution as ``d``.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> params(d)
(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)

```

"""
params(d::SGS) = (d.E, d.b, d.g)

# Required methods for a UnivariateDistribution type
# See https://juliastats.github.io/Distributions.jl/latest/extends.html#Univariate-Distribution-1
# rand(d::Power)
# sampler(d::Power)
# pdf(d::Power, x::Real)
# logpdf(d::Power, x::Real)
# cdf(d::Power, x::Real)
# quantile(d::Power, q::Real)
# minimum(d::Power)
# maximum(d::Power)
# insupport(d::Power, x::Real)


# Sampling

"""

    invitp(d::SGS)

Use a linear interpolation to produce an inversely transformed sample of the
SGS distribution.

# External links

* [Inverse transform sampling on Wikipedia](https://en.wikipedia.org/wiki/Inverse_transform_sampling)

"""
function invitp(d::SGS)
    xarr = collect(-10:0.01:10)
    A = cdf.(d, xarr)
    itp = LinearInterpolation(A, xarr, extrapolation_bc=Line())
    return itp
end


"""

    rand(d::SGS)

Generate a scalar sample from the SGS distribution `d`.

"""
function rand(d::SGS)
    return rand(GLOBAL_RNG, d)
end

"""

    rand(d::SGS, n::Int64)

Generate an array of `n` samples from the SGS distribution `d`.

"""
function rand(d::SGS, n::Int64)
    return rand(GLOBAL_RNG, d, n)
end

"""

    rand(d::SGS)

Generate a scalar sample from the SGS distribution `d`.

"""
function rand(rng::AbstractRNG, d::SGS)
    itp = invitp(d)
    return itp(rand(rng))
end

"""

    rand(rng::AbstractRNG, d::SGS, n::Int64)

Generate an array of `n` samples from the SGS distribution `d`.

"""
function rand(rng::AbstractRNG, d::SGS, n::Int64)
    itp = invitp(d)
    return itp(rand(rng, n))
end


# Statistics


"""

    mean(d::SGS)

Compute the mean of the SGS distribution.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> mean(d)
0.0

```

"""
function mean(d::SGS)
    return 0.0
end

"""

    var(d::SGS)

Compute the variance of the SGS distribution.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> var(d)
1.0

```

"""
function var(d::SGS)
    (E, b, g) = params(d)
    return (g^2 + b^2) / (2 - E^2)
end


"""

    std(d::SGS)

Compute the standard deviation of the SGS distribution.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> std(d)
1.0

```

"""
function std(d::SGS)
    return sqrt(var(d))
end


"""

    skewness(d::SGS)

Compute the skewness of the SGS distribution.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> skewness(d)
1.0

```

"""
function skewness(d::SGS)
    (E, b, g) = params(d)
    return moment_constraint(3, E) ? ((2.0 * E) / (1.0 - E^2)) * (g / sqrt(var(d))) : NaN
end


"""

    kurtosis(d::SGS; correction::Bool=true)

Compute the excess kurtosis of the SGS distribution. Excess kurtosis is returned by default with `correction=true`.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> kurtosis(d, true)
4.999999999999998

```

"""
function kurtosis(d::SGS; excess::Bool=true)
    (E, b, g) = params(d)
    if moment_constraint(4, E)
        term1 = 1.5 * ((1.0 - E^2) / (1.0 - 1.5 * E^2)) * skewness(d)^2
        term2 = (3.0 * E^2) / (1.0 - 1.5*E^2)
        kurt = term1 + term2
        if excess
            return kurt         # excess kurtosis
        else
            return kurt + 3     # kurtosis
        end
    else
        return NaN  # kurtosis does not exist
    end
    # Skewness-kurtosis constraint.
    # if excess_kurtosis ≤ 1.5 * skewness(d)^2
    #     return NaN
    # else
    #     return excess_kurtosis
end


# Evaluation

"""

    pdf(d::SGS, x::Real; norm::Bool=true)

Compute the pdf of the SGS distribution. The SGS pdf is normalized by default,
i.e., `norm=true`. If `norm` is set to `false`, then the SGS normalization
constant will not be applied.

Here, the SGS probability density function (pdf) is given in the form of [Sardeshmukh et al. (2015)](https://doi.org/10.1175/JCLI-D-15-0020.1) as

```math
p(x) = \\frac{1}{\\mathcal{N}} \\left[ \\left( Ex + g \\right)^{2} + b^{2}  \\right]^{-\\left(1 + \\left(1/E^{2}\\right)\\right)} \\mathrm{exp} \\! \\left[ \\frac{2g}{E^{2}b} \\mathrm{arctan} \\! \\left(\\frac{Ex + g}{b} \\right) \\right]
```

where the normalization constant, ``\\mathcal{N}``, is

```math
\\mathcal{N} = \\frac{2 \\pi \\nu^{1/2} \\left( 2b \\right)^{-\\left(2 \\nu + 1 \\right)} \\Gamma \\! \\left(2\\nu + 1\\right)}{\\Gamma \\! \\left( \\nu + 1 - iq/2 \\right) \\Gamma \\! \\left( \\nu + 1 + iq/2 \\right)}
```

See also:

>Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. *Journal of Climate*, **28**, 9166–9187, [https://doi.org/10.1175/JCLI-D-15-0020.1](https://doi.org/10.1175/JCLI-D-15-0020.1).

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> pdf(d, 0.0)
0.46210984589674214

```

"""
function pdf(d::SGS, x::Real; norm::Bool=true)
    return _pdf(d, x, norm)
end

"""

    pdf(d::SGS, x::AbstractArray; norm::Bool=true)

Compute the pdf of the SGS distribution. The SGS pdf is normalized by default,
i.e., `norm=true`. If `norm` is set to `false`, then the SGS normalization
constant will not be applied.

Here, the SGS probability density function (pdf) is given in the form of [Sardeshmukh et al. (2015)](https://doi.org/10.1175/JCLI-D-15-0020.1) as

```math
p(x) = \\frac{1}{\\mathcal{N}} \\left[ \\left( Ex + g \\right)^{2} + b^{2}  \\right]^{-\\left(1 + \\left(1/E^{2}\\right)\\right)} \\mathrm{exp} \\! \\left[ \\frac{2g}{E^{2}b} \\mathrm{arctan} \\! \\left(\\frac{Ex + g}{b} \\right) \\right]
```

where the normalization constant, ``\\mathcal{N}``, is

```math
\\mathcal{N} = \\frac{2 \\pi \\nu^{1/2} \\left( 2b \\right)^{-\\left(2 \\nu + 1 \\right)} \\Gamma \\! \\left(2\\nu + 1\\right)}{\\Gamma \\! \\left( \\nu + 1 - iq/2 \\right) \\Gamma \\! \\left( \\nu + 1 + iq/2 \\right)}
```

See also:

>Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. *Journal of Climate*, **28**, 9166–9187, [https://doi.org/10.1175/JCLI-D-15-0020.1](https://doi.org/10.1175/JCLI-D-15-0020.1).

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = [-0.2, 0.4, 1.5, 5.0]
4-element Array{Float64,1}:
 -0.2
  0.4
  1.5
  5.0

julia> pdf(d, x)
4-element Array{Float64,1}:
 0.4821612875311956
 0.3552310712921366
 0.091209594116695
 0.0011832665173870727


```

"""
function pdf(d::SGS, x::AbstractArray; norm::Bool=true)
    pdfs = Array{Float64}(undef, length(x))
    for (indx, num) in enumerate(x)
        pdfs[indx] = _pdf(d, num, norm)
    end
    return pdfs
end


function pdf(E::Real, b::Real, g::Real, x::Real)
    return _pdf(SGS(E, b, g), x, false)
end


function _pdf(d::SGS, x::Real, norm::Bool)
    (E, b, g) = params(d)

    @debug "E = $E"
    @debug "b = $b"
    @debug "g = $g"

    # if E < 0.05
    #     return _pdf_E_lim(d, x)
    # end

    # Calculate the log of the normalization constant.
    if norm
        log_N = logN(E, b, g)
    else
        log_N = 0.0
    end

    # Calculate the log pdf.
    q = (2.0 * g) / (E^2 * b)
    @debug "p(x) q = $q"
    atan_term = atan((E * x + g) / b)
    @debug "p(x) atan_term = $atan_term"
    log_px = -log_N + q * atan_term - (1.0 + 1.0/E^2) * log((E * x + g)^2 + b^2)
    @debug "log(p(x)) = $log_px"

    # Calculate the pdf.
    px = exp(log_px)
    @debug "p(x) = $px"

    return px

end


function _pdf(E::Real, b::Real, g::Real, x::Real)

    @debug "E = $E"
    @debug "b = $b"
    @debug "g = $g"

    # Calculate the log of the normalization constant.
    # if norm
    #     log_N = logN(E, b, g)
    # else
    #     log_N = 0.0
    # end

    log_N = 0.0

    # Calculate the log pdf.
    q = (2.0 * g) / (E^2 * b)
    @debug "p(x) q = $q"
    atan_term = atan((E * x + g) / b)
    @debug "p(x) atan_term = $atan_term"
    log_px = -log_N + q * atan_term - (1.0 + 1.0/E^2) * log((E * x + g)^2 + b^2)
    @debug "log(p(x)) = $log_px"

    # Calculate the pdf.
    px = exp(log_px)
    @debug "p(x) = $px"

    return px

end


function _pdf_E_lim(d::SGS, x::Real)
    (E, b, g) = params(d)
    px = 1.0 / sqrt(2*π) * exp((-x^2 / 2.0) + ((E * g * x) / 3.0) * (x^2 - 3.0))
    return px
end

"""

    pdf_xmax(d::SGS)

Compute the location of the unique pdf maximum.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> pdf_xmax(d)
-0.22000000000000003

```

"""
function pdf_xmax(d::SGS)
    (E, b, g) = params(d)
    return -1.0 * E * g / (1.0 + E^2)
end


"""

    logN(E::Real, b::Real, g::Real)

Compute the log of the normalization constant `N` of the SGS distribution.

# Examples
```jldoctest
julia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> logN(d.E, d.b, d.g)
-0.07844022339848955

```

"""
function logN(E::Real, b::Real, g::Real)

    # Calculate parameters for the normalization constant.
    ν = 1.0 / E^2
    q = (2.0 * g * ν) / b
    @debug "ν = $ν"
    @debug "q = $q"

    # Calculate first two terms in numerator, taking logs for numerical stability.
    term1 = log(2.0 * π * sqrt(ν))
    term2 = log((2.0 * b)^(-(2.0 * ν + 1.0)))
    @debug "log(2.0 * π * sqrt(ν)) = $term1"
    @debug "log((2.0 * b)^(1.0 - 2.0 * ν)) = $term2"

    # Calculate parameters for the gamma function.
    logΓ0 = loggamma(2.0 * ν + 1.0)
    logΓ1 = loggamma(ν + 1.0 - (1im * q / 2.0))
    logΓ2 = loggamma(ν + 1.0 + (1im * q / 2.0))
    @debug "log(Γ0) = $logΓ0"
    @debug "log(Γ1) = $logΓ1"
    @debug "log(Γ2) = $logΓ2"

    # Find the log of the normalization constant.
    logN = real(term1 + term2 + logΓ0 - logΓ1 - logΓ2)
    @debug "log(N) = $logN"
    # N = exp(logN)
    # @debug "N = $N"
    #
    # @debug "Returning logN..."
    return logN

end


"""

    logpdf(d::SGS, x::AbstractArray)

Calculate the log of the probability density function of the SGS distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
8-element Array{Float64,1}:
 -10.0
  -1.0
  -0.5
   0.0
   0.5
   1.0
   2.0
   5.0

julia> logpdf(d, x)
8-element Array{Float64,1}:
 -15.504184293152903
  -1.3393373142985014
  -0.8053427806087643
  -0.7719526544799672
  -1.1296521694192523
  -1.7130802377432286
  -3.099847429101194
  -6.739476429937023

```

"""
function logpdf(d::SGS, x::AbstractArray)
    lpdfs = Array{Float64}(undef, length(x))
    for (indx, num) in enumerate(x)
        lpdfs[indx] = log(pdf(d, num))
    end
    return lpdfs
end


"""

    logpdf(d::SGS, x::Real)

Calculate the log of the cumulative distribution function of the SGS
distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = 0.5
0.5

julia> logpdf(d, x)
-1.1296521694192523

```

"""
function logpdf(d::SGS, x::Real)
    return log(pdf(d, x))
end


"""

    cdf(d::SGS, x::AbstractArray)

Calculate the cumulative distribution function of the SGS distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
8-element Array{Float64,1}:
 -10.0
  -1.0
  -0.5
   0.0
   0.5
   1.0
   2.0
   5.0

julia> cdf(d, x)
8-element Array{Float64,1}:
 2.7059681332492457e-7
 0.1268532553530455
 0.3074639317904149
 0.5437968749444915
 0.7430905390656029
 0.86706141301943
 0.9654409769867174
 0.9986564347567426

```

"""
function cdf(d::SGS, x::AbstractArray)
    cdfs = Array{Float64}(undef, length(x))
    for (indx, num) in enumerate(x)
        cdfs[indx] = cdf(d, num)
    end
    return cdfs
end


"""

    cdf(d::SGS, x::Real)

Calculate the cumulative distribution function of the SGS distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = 0.5
0.5

julia> cdf(d, x)
0.7430905390656029

```

"""
function cdf(d::SGS, x::Real)
    integral, error = quadgk(x -> pdf(d, x), -Inf, x)
    return integral
end


# function cdf(d::SGS, x::Real)
#
#     dx = 0.01
#     σval = collect(-10:dx:x)
#
#     cdf = 0.0
#     for i in σval
#         cdf += pdf(d, i) * dx
#     end
#
#     return cdf
#
# end
#
# function cdf(d::SGS, x::AbstractArray)
#
#     dx = 0.01
#     cdf = zeros(length(x))
#     for i=1:length(x)
#         for xval in -10:dx:x[i]
#             cdf[i] += pdf(d, xval) * dx
#         end
#     end
#     return cdf
#
# end


"""

    logcdf(d::SGS, x::AbstractArray)

Calculate the log of the cumulative distribution function of the SGS
distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
8-element Array{Float64,1}:
 -10.0
  -1.0
  -0.5
   0.0
   0.5
   1.0
   2.0
   5.0

julia> logcdf(d, x)
8-element Array{Float64,1}:
 -15.12263589760972
  -2.064724330254346
  -1.1793974936056804
  -0.6091794935003683
  -0.2969373856106992
  -0.14264547077775455
  -0.035170311051895865
  -0.0013444686363078492

```

"""
function logcdf(d::SGS, x::AbstractArray)
    lcdfs = Array{Float64}(undef, length(x))
    for (indx, num) in enumerate(x)
        lcdfs[indx] = log(cdf(d, num))
    end
    return lcdfs
end


"""

    logcdf(d::SGS, x::Real)

Calculate the log of the cumulative distribution function of the SGS
distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = 0.5
0.5

julia> logcdf(d, x)
-0.2969373856106992

```

"""
function logcdf(d::SGS, x::Real)
    integral, error = quadgk(x -> pdf(d, x), -Inf, x)
    return log(integral)
end

# function logcdf(d::SGS, x::Real)
#     return log(cdf(d, x))
# end


"""

    ccdf(d::SGS, x::AbstractArray)

Calculate the complementary cumulative distribution function of the SGS
distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
8-element Array{Float64,1}:
 -10.0
  -1.0
  -0.5
   0.0
   0.5
   1.0
   2.0
   5.0

julia> ccdf(d, x)
8-element Array{Float64,1}:
 0.9999997294031867
 0.8731467446469545
 0.6925360682095851
 0.45620312505550853
 0.25690946093439715
 0.13293858698057004
 0.03455902301328262
 0.0013435652432574052

```

"""
function ccdf(d::SGS, x::AbstractArray)
    ccdfs = Array{Float64}(undef, length(x))
    for (indx, num) in enumerate(x)
        ccdfs[indx] = 1.0 - cdf(d, num)
    end
    return ccdfs
end

# function ccdf(d::SGS, x::AbstractArray)
#     ccdfs = Array{Float64}(undef, length(x))
#     for (indx, num) in enumerate(x)
#         ccdfs[indx] = ccdf(d, num)
#     end
#     return ccdfs
# end


"""

    ccdf(d::SGS, x::Real)

Calculate the complementary cumulative distribution function of the SGS
distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = 0.5
0.5

julia> ccdf(d, x)
0.25690946093439715

```

"""
function ccdf(d::SGS, x::Real)
    integral, error = quadgk(x -> pdf(d, x), -Inf, x)
    return 1.0 - integral
end

# function ccdf(d::SGS, x::Real)
#     return 1.0 - cdf(d, x)
# end


"""

    logccdf(d::SGS, x::AbstractArray)

Calculate the log of the complementary cumulative distribution function of the
SGS distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]
8-element Array{Float64,1}:
 -10.0
  -1.0
  -0.5
   0.0
   0.5
   1.0
   2.0
   5.0

julia> logccdf(d, x)
8-element Array{Float64,1}:
 -2.70596849901216e-7
 -0.1356516448893756
 -0.3673949582198083
 -0.7848171189678753
 -1.3590315482404094
 -2.017868009426468
 -3.365086604737457
 -6.612428568931222

```

"""
function logccdf(d::SGS, x::AbstractArray)
    ccdfs = Array{Float64}(undef, length(x))
    for (indx, num) in enumerate(x)
        ccdfs[indx] = log(1.0 - cdf(d, num))
    end
    return ccdfs
end

# function logccdf(d::SGS, x::AbstractArray)
#     logccdfs = Array{Float64}(undef, length(x))
#     for (indx, num) in enumerate(x)
#         logccdfs[indx] = log(ccdf(d, num))
#     end
#     return logccdfs
# end


"""

    logccdf(d::SGS, x::Real)

Calculate the log of the complementary cumulative distribution function of the
SGS distribution.

# Examples
```jldoctest
julia> d = fit_mm(SGS, 1, 1, 5)
SGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)

julia> x = 0.5
0.5

julia> logccdf(d, x)
-1.3590315482404094

```

"""
function logccdf(d::SGS, x::Real)
    integral, error = quadgk(x -> pdf(d, x), -Inf, x)
    return log(1.0 - integral)
end

# function logccdf(d::SGS, x::Real)
#     return log(ccdf(d, x))
# end
