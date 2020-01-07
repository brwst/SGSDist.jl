__precompile__()

module SGSDist

    using Logging
    import DelimitedFiles: readdlm
    import Random: rand, randn, AbstractRNG, GLOBAL_RNG

    import DifferentialEquations: SDEProblem, solve, EulerHeun, EM, RKMil
    import Distributions: @distr_support, @check_args, ContinuousUnivariateDistribution, Normal, params, pdf, cdf, ccdf, logpdf, logcdf, logccdf, fit, mean, var, skewness, kurtosis
    import Interpolations: LinearInterpolation, Line
    import QuadGK: quadgk
    import SpecialFunctions: loggamma
    import StatsBase: var, skewness, kurtosis

    export SGS, params, mean, var, std, skewness, kurtosis, pdf, pdf_xmax, logN, logpdf, cdf, logcdf, ccdf, logccdf, fit, CAM1D, Hasselmann1D, RK4, Euler, fdt, rand, air, epsilon_condition, epsilon_constraint, ks_constraint, small_E_constraint, fit_mm, param_E, param_b, param_g, invitp

    include("constraints.jl")
    include("continuous.jl")
    include("fit.jl")
    include("markov.jl")
    include("sampledata.jl")

    # precompile hints
    # d = SGS()
    # rand(d, 1)

end
