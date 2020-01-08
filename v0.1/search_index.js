var documenterSearchIndex = {"docs":
[{"location":"api/#Index-1","page":"Index","title":"Index","text":"","category":"section"},{"location":"api/#","page":"Index","title":"Index","text":"","category":"page"},{"location":"library/#Library-1","page":"Library","title":"Library","text":"","category":"section"},{"location":"library/#","page":"Library","title":"Library","text":"CurrentModule = SGSDist","category":"page"},{"location":"library/#","page":"Library","title":"Library","text":"SGSDist.SGS\nSGSDist.params\nSGSDist.mean\nSGSDist.var\nSGSDist.std\nSGSDist.skewness\nSGSDist.kurtosis\nSGSDist.pdf\nSGSDist.logpdf\nSGSDist.cdf\nSGSDist.logcdf\nSGSDist.ccdf\nSGSDist.logccdf\nSGSDist.fit\nSGSDist.CAM1D\nSGSDist.Hasselmann1D\nSGSDist.fdt\nSGSDist.rand\nSGSDist.air","category":"page"},{"location":"library/#SGSDist.SGS","page":"Library","title":"SGSDist.SGS","text":"SGS(E, b, g)\n\nThe stochastically generated skewed (SGS) distribution with parameters E, b and g, derived by Sardeshmukh and Sura (2009) and following the form of Sardeshmukh et al. (2015). The SGS probability density function (pdf) is\n\np(x) = frac1mathcalN left left( Ex + g right)^2 + b^2  right^-left(1 + left(1E^2right)right) mathrmexp  left frac2gE^2b mathrmarctan  left(fracEx + gb right) right\n\nwhere the normalization constant, mathcalN, is given as\n\nmathcalN = frac2 pi nu^12 left( 2b right)^-left(2 nu + 1 right) Gamma  left(2nu + 1right)Gamma  left( nu + 1 - iq2 right) Gamma  left( nu + 1 + iq2 right)\n\nSee also:\n\nSardeshmukh, P. D., and P. Sura, 2009: Reconciling Non-Gaussian Climate Statistics with Linear Dynamics. Journal of Climate, 22, 1193–1207, https://doi.org/10.1175/2008JCLI2358.1.Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. Journal of Climate, 28, 9166–9187, https://doi.org/10.1175/JCLI-D-15-0020.1.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\n\n\n\n\n\n","category":"type"},{"location":"library/#StatsBase.params","page":"Library","title":"StatsBase.params","text":"params(d::SGS)\n\nReturn a tuple of SGS parameters E, b and g. Let d be a distribution of type D, then D(params(d)...) will construct exactly the same distribution as d.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> params(d)\n(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Statistics.mean","page":"Library","title":"Statistics.mean","text":"mean(d::SGS)\n\nCompute the mean of the SGS distribution.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> mean(d)\n0.0\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Statistics.var","page":"Library","title":"Statistics.var","text":"var(d::SGS)\n\nCompute the variance of the SGS distribution.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> var(d)\n1.0\n\n\n\n\n\n\n","category":"function"},{"location":"library/#SGSDist.std","page":"Library","title":"SGSDist.std","text":"std(d::SGS)\n\nCompute the standard deviation of the SGS distribution.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> std(d)\n1.0\n\n\n\n\n\n\n","category":"function"},{"location":"library/#StatsBase.skewness","page":"Library","title":"StatsBase.skewness","text":"skewness(d::SGS)\n\nCompute the skewness of the SGS distribution.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> skewness(d)\n1.0\n\n\n\n\n\n\n","category":"function"},{"location":"library/#StatsBase.kurtosis","page":"Library","title":"StatsBase.kurtosis","text":"kurtosis(d::SGS; correction::Bool=true)\n\nCompute the excess kurtosis of the SGS distribution. Excess kurtosis is returned by default with correction=true.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> kurtosis(d, true)\n4.999999999999998\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Distributions.pdf","page":"Library","title":"Distributions.pdf","text":"pdf(d::SGS, x::Real; norm::Bool=true)\n\nCompute the pdf of the SGS distribution. The SGS pdf is normalized by default, i.e., norm=true. If norm is set to false, then the SGS normalization constant will not be applied.\n\nHere, the SGS probability density function (pdf) is given in the form of Sardeshmukh et al. (2015) as\n\np(x) = frac1mathcalN left left( Ex + g right)^2 + b^2  right^-left(1 + left(1E^2right)right) mathrmexp  left frac2gE^2b mathrmarctan  left(fracEx + gb right) right\n\nwhere the normalization constant, mathcalN, is\n\nmathcalN = frac2 pi nu^12 left( 2b right)^-left(2 nu + 1 right) Gamma  left(2nu + 1right)Gamma  left( nu + 1 - iq2 right) Gamma  left( nu + 1 + iq2 right)\n\nSee also:\n\nSardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. Journal of Climate, 28, 9166–9187, https://doi.org/10.1175/JCLI-D-15-0020.1.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> pdf(d, 0.0)\n0.46210984589674214\n\n\n\n\n\n\npdf(d::SGS, x::AbstractArray; norm::Bool=true)\n\nCompute the pdf of the SGS distribution. The SGS pdf is normalized by default, i.e., norm=true. If norm is set to false, then the SGS normalization constant will not be applied.\n\nHere, the SGS probability density function (pdf) is given in the form of Sardeshmukh et al. (2015) as\n\np(x) = frac1mathcalN left left( Ex + g right)^2 + b^2  right^-left(1 + left(1E^2right)right) mathrmexp  left frac2gE^2b mathrmarctan  left(fracEx + gb right) right\n\nwhere the normalization constant, mathcalN, is\n\nmathcalN = frac2 pi nu^12 left( 2b right)^-left(2 nu + 1 right) Gamma  left(2nu + 1right)Gamma  left( nu + 1 - iq2 right) Gamma  left( nu + 1 + iq2 right)\n\nSee also:\n\nSardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. Journal of Climate, 28, 9166–9187, https://doi.org/10.1175/JCLI-D-15-0020.1.\n\nExamples\n\njulia> d = SGS(0.6236095644623235, 1.1709106481844573, 0.48997894350611143)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = [-0.2, 0.4, 1.5, 5.0]\n4-element Array{Float64,1}:\n -0.2\n  0.4\n  1.5\n  5.0\n\njulia> pdf(d, x)\n4-element Array{Float64,1}:\n 0.4821612875311956\n 0.3552310712921366\n 0.091209594116695\n 0.0011832665173870727\n\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Distributions.logpdf","page":"Library","title":"Distributions.logpdf","text":"logpdf(d::SGS, x::AbstractArray)\n\nCalculate the log of the probability density function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]\n8-element Array{Float64,1}:\n -10.0\n  -1.0\n  -0.5\n   0.0\n   0.5\n   1.0\n   2.0\n   5.0\n\njulia> logpdf(d, x)\n8-element Array{Float64,1}:\n -15.504184293152903\n  -1.3393373142985014\n  -0.8053427806087643\n  -0.7719526544799672\n  -1.1296521694192523\n  -1.7130802377432286\n  -3.099847429101194\n  -6.739476429937023\n\n\n\n\n\n\nlogpdf(d::SGS, x::Real)\n\nCalculate the log of the cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = 0.5\n0.5\n\njulia> logpdf(d, x)\n-1.1296521694192523\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Distributions.cdf","page":"Library","title":"Distributions.cdf","text":"cdf(d::SGS, x::AbstractArray)\n\nCalculate the cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]\n8-element Array{Float64,1}:\n -10.0\n  -1.0\n  -0.5\n   0.0\n   0.5\n   1.0\n   2.0\n   5.0\n\njulia> cdf(d, x)\n8-element Array{Float64,1}:\n 2.7059681332492457e-7\n 0.1268532553530455\n 0.3074639317904149\n 0.5437968749444915\n 0.7430905390656029\n 0.86706141301943\n 0.9654409769867174\n 0.9986564347567426\n\n\n\n\n\n\ncdf(d::SGS, x::Real)\n\nCalculate the cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = 0.5\n0.5\n\njulia> cdf(d, x)\n0.7430905390656029\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Distributions.logcdf","page":"Library","title":"Distributions.logcdf","text":"logcdf(d::SGS, x::AbstractArray)\n\nCalculate the log of the cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]\n8-element Array{Float64,1}:\n -10.0\n  -1.0\n  -0.5\n   0.0\n   0.5\n   1.0\n   2.0\n   5.0\n\njulia> logcdf(d, x)\n8-element Array{Float64,1}:\n -15.12263589760972\n  -2.064724330254346\n  -1.1793974936056804\n  -0.6091794935003683\n  -0.2969373856106992\n  -0.14264547077775455\n  -0.035170311051895865\n  -0.0013444686363078492\n\n\n\n\n\n\nlogcdf(d::SGS, x::Real)\n\nCalculate the log of the cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = 0.5\n0.5\n\njulia> logcdf(d, x)\n-0.2969373856106992\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Distributions.ccdf","page":"Library","title":"Distributions.ccdf","text":"ccdf(d::SGS, x::AbstractArray)\n\nCalculate the complementary cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]\n8-element Array{Float64,1}:\n -10.0\n  -1.0\n  -0.5\n   0.0\n   0.5\n   1.0\n   2.0\n   5.0\n\njulia> ccdf(d, x)\n8-element Array{Float64,1}:\n 0.9999997294031867\n 0.8731467446469545\n 0.6925360682095851\n 0.45620312505550853\n 0.25690946093439715\n 0.13293858698057004\n 0.03455902301328262\n 0.0013435652432574052\n\n\n\n\n\n\nccdf(d::SGS, x::Real)\n\nCalculate the complementary cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = 0.5\n0.5\n\njulia> ccdf(d, x)\n0.25690946093439715\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Distributions.logccdf","page":"Library","title":"Distributions.logccdf","text":"logccdf(d::SGS, x::AbstractArray)\n\nCalculate the log of the complementary cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = [-10.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 5.0]\n8-element Array{Float64,1}:\n -10.0\n  -1.0\n  -0.5\n   0.0\n   0.5\n   1.0\n   2.0\n   5.0\n\njulia> logccdf(d, x)\n8-element Array{Float64,1}:\n -2.70596849901216e-7\n -0.1356516448893756\n -0.3673949582198083\n -0.7848171189678753\n -1.3590315482404094\n -2.017868009426468\n -3.365086604737457\n -6.612428568931222\n\n\n\n\n\n\nlogccdf(d::SGS, x::Real)\n\nCalculate the log of the complementary cumulative distribution function of the SGS distribution.\n\nExamples\n\njulia> d = fit_mm(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = 0.5\n0.5\n\njulia> logccdf(d, x)\n-1.3590315482404094\n\n\n\n\n\n\n","category":"function"},{"location":"library/#StatsBase.fit","page":"Library","title":"StatsBase.fit","text":"fit(d::Type{SGS}, x::AbstractArray)\n\nFit an SGS distribution to the time series, x.\n\nExamples\n\njulia> fit(SGS, air())\nSGS{Float64}(E=0.36623331624359673, b=0.7133796999284913, g=-1.1648873051087967)\n\n\n\n\n\n\nfit(d::Type{SGS}, variance::Real, skew::Real, kurt::Real)\n\nFit an SGS distribution to the desired variance, skewness and kurtosis.\n\nExamples\n\njulia> fit(SGS, 1.0, 1.0, 5.0)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\n\n\n\n\n\n","category":"function"},{"location":"library/#SGSDist.CAM1D","page":"Library","title":"SGSDist.CAM1D","text":"CAM1D(d::SGS, n::Integer; dt::Real=1/24, λ::Real=1.0)\n\nCreate a timeseries using the one-dimensional CAM noise model of Sardeshmukh and Sura (2009) for a particular SGS distribution.\n\nSee also:\n\nSardeshmukh, P. D., and P. Sura, 2009: Reconciling Non-Gaussian Climate Statistics with Linear Dynamics. Journal of Climate, 22, 1193–1207, https://doi.org/10.1175/2008JCLI2358.1.Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. Journal of Climate, 28, 9166–9187, https://doi.org/10.1175/JCLI-D-15-0020.1.\n\nExamples\n\njulia> d = fit(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)\n\njulia> x = CAM1D(d, 1000, seed=42)\nretcode: Success\nInterpolation: 1st order linear\nt: 1001-element Array{Float64,1}:\n  0.0\n  0.041666666666666664\n  0.08333333333333333\n  0.125\n  0.16666666666666666\n  0.20833333333333331\n  0.24999999999999997\n  0.29166666666666663\n  0.3333333333333333\n  0.375\n  ⋮\n 41.333333333333165\n 41.37499999999983\n 41.416666666666494\n 41.45833333333316\n 41.49999999999982\n 41.54166666666649\n 41.58333333333315\n 41.624999999999815\n 41.625\nu: 1001-element Array{Array{Float64,1},1}:\n [0.0]\n [-0.05538507236046931]\n [0.3759078513765659]\n [0.31482891888314446]\n [0.2586995978771358]\n [0.0063144005340343146]\n [0.0819076515551697]\n [0.07880240481554046]\n [-0.09903540075387843]\n [-0.0867782526753085]\n ⋮\n [-0.048583800782025996]\n [0.30770620643597935]\n [0.6914165604620269]\n [1.2048828846162343]\n [0.6995798768325538]\n [0.5417553247799438]\n [0.593919650292867]\n [0.6785936302181321]\n [0.6785940957553152]\n\n\n\n\n\n\n","category":"function"},{"location":"library/#SGSDist.Hasselmann1D","page":"Library","title":"SGSDist.Hasselmann1D","text":"Hasselmann1D(n::Integer; dt::Real=1/24, λ::Real=1.0, seed::Real=-1)\n\nCreate a timeseries using the one-dimensional noise model of Hasselmann (1976).\n\nSee also:\n\nHasselmann, K., 1976: Stochastic climate models Part I. Theory. Tellus, 28, 473–485, https://doi.org/10.1111/j.2153-3490.1976.tb00696.x.\n\nExamples\n\njulia> Hasselmann1D(1000, seed=42)\nretcode: Success\nInterpolation: 1st order linear\nt: 1001-element Array{Float64,1}:\n  0.0\n  0.041666666666666664\n  0.08333333333333333\n  0.125\n  0.16666666666666666\n  0.20833333333333331\n  0.24999999999999997\n  0.29166666666666663\n  0.3333333333333333\n  0.375\n  ⋮\n 41.333333333333165\n 41.37499999999983\n 41.416666666666494\n 41.45833333333316\n 41.49999999999982\n 41.54166666666649\n 41.58333333333315\n 41.624999999999815\n 41.625\nu: 1001-element Array{Array{Float64,1},1}:\n [0.0]\n [-0.13005535542166108]\n [0.04828669603525988]\n [0.658691166050971]\n [0.4766431404349038]\n [0.2527905702593073]\n [0.5063085613562422]\n [0.4424101651546515]\n [0.42457583766050366]\n [0.17019337828909467]\n ⋮\n [0.04255470162528974]\n [-0.009396010196968833]\n [-0.4700188956712661]\n [-0.6640224731670598]\n [-0.6412697125734018]\n [-0.8226707386234939]\n [-0.9240097641730844]\n [-0.3550868399241892]\n [-0.35508653697280324]\n\n\n\n\n\n\n","category":"function"},{"location":"library/#SGSDist.fdt","page":"Library","title":"SGSDist.fdt","text":"fdt(λ::Real, ν::Real, dt::Real)\n\nCompute the variance that satisfies the fluctuation-dissipation theorem.\n\nExternal links\n\nFluctuation-dissipation theorem on Wikipedia\n\nExamples\n\njulia> fdt(1.0, 1.0, 0.1)\n1.378404875209022\n\n\n\n\n\n\n","category":"function"},{"location":"library/#Base.rand","page":"Library","title":"Base.rand","text":"rand(d::SGS)\n\nGenerate a scalar sample from the SGS distribution d.\n\n\n\n\n\nrand(d::SGS, n::Int64)\n\nGenerate an array of n samples from the SGS distribution d.\n\n\n\n\n\nrand(d::SGS)\n\nGenerate a scalar sample from the SGS distribution d.\n\n\n\n\n\nrand(rng::AbstractRNG, d::SGS, n::Int64)\n\nGenerate an array of n samples from the SGS distribution d.\n\n\n\n\n\n","category":"function"},{"location":"library/#SGSDist.air","page":"Library","title":"SGSDist.air","text":"air()\n\nProduce an array of sample climate data for testing. This sample time series is taken from the NCEP/NCAR Reanalysis 1 dataset, representing the 925 hPa air temperature standardized anomalies of all  DJF winters from 1948-2018 at 50.0°N, 235.0°E.\n\nSee also:\n\nKalnay, E., and Coauthors, 1996: The NCEP/NCAR 40-year reanalysis project. Bulletin of the American Meteorological Society, 77, 437–471, https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2.\n\nExamples\n\njulia> air()\n6408-element Array{Float64,1}:\n  0.67547214\n  0.8308915\n  0.19012427\n -0.20711477\n -0.10181273\n  0.7421673\n  0.22026125\n  0.050515447\n  0.043254483\n  0.044615462\n  ⋮\n  0.7591071\n  0.54696834\n  0.3421846\n  0.044188347\n  0.2433472\n  0.4460243\n  1.4941149\n  0.49388295\n -0.052520715\n\n\n\n\n\n\n","category":"function"},{"location":"examples/#Examples-1","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Create-an-SGS-distribution-1","page":"Examples","title":"Create an SGS distribution","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"One may create an SGS distribution by fitting the distribution to known statistics by the method of moments. For example, to create an SGS distribution with mean zero, a variance and skewness of 1 and a kurtosis of 5, issue:","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> using SGSDist\n\njulia> d = fit(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Here, d is an SGS distribution that has been fitted to the moments specified with parameters E, g and b.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"By loading the Plots.jl package, one can plot this SGS pdf:","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> using Plots\n\njulia> x = collect(-6:0.01:6)\n1201-element Array{Float64,1}:\n -6.0\n -5.99\n -5.98\n -5.97\n -5.96\n -5.95\n -5.94\n  ⋮   \n  5.94\n  5.95\n  5.96\n  5.97\n  5.98\n  5.99\n  6.0\n\njulia> plot(x, pdf.(d, x))\n","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: Sample SGS pdf)","category":"page"},{"location":"examples/#Fit-an-SGS-distribution-from-data-1","page":"Examples","title":"Fit an SGS distribution from data","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"An SGS distribution may be also be fit from climate data. For example, one may create an SGS distribution by fitting the sample reanalysis time series provided by the function air():","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> air()\n6408-element Array{Float64,1}:\n  0.67547214\n  0.8308915  \n  0.19012427\n -0.20711477\n -0.10181273\n  0.7421673  \n  0.22026125\n  ⋮          \n  0.3421846  \n  0.044188347\n  0.2433472  \n  0.4460243  \n  1.4941149  \n  0.49388295\n -0.052520715\n\njulia> d = fit(SGS, air())\nSGS{Float64}(E=0.36623331624359673, b=0.7133796999284913, g=-1.1648873051087967)\n\njulia> plot(x, pdf.(d, x))\n\njulia> using StatsPlots\n\njulia> hist!(x)\n","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: Sample air() pdf)","category":"page"},{"location":"examples/#Markov-processes-1","page":"Examples","title":"Markov processes","text":"","category":"section"},{"location":"examples/#Correlated-and-additive-multiplicative-(CAM)-noise-model-1","page":"Examples","title":"Correlated and additive multiplicative (CAM) noise model","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"One may create a time series via Markov process that represents the statistics of a given SGS distribution. First, specify the SGS distribution – here we create an SGS distribution with variance and skewness of 1 and (excess) kurtosis of 5.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> d = fit(SGS, 1, 1, 5)\nSGS{Float64}(E=0.6236095644623235, b=1.1709106481844573, g=0.48997894350611143)","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"We can produce a time series using the CAM1D() function, which corresponds to the one-dimensional correlated additive and multiplicative (CAM) noise model described in Sardeshmukh and Sura (2009). The CAM noise model is able to reproduce non-Gaussian statistics that are seen in climate observations of daily atmospheric variables.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> n = 1000  # length of time series\n1000\n\njulia> sol = CAM1D(d, n)\nretcode: Success\nInterpolation: 1st order linear\nt: 1001-element Array{Float64,1}:\n  0.0                 \n  0.041666666666666664\n  0.08333333333333333\n  0.125               \n  0.16666666666666666\n  ⋮                   \n 41.49999999999982    \n 41.54166666666649    \n 41.58333333333315    \n 41.624999999999815   \n 41.625               \nu: 1001-element Array{Array{Float64,1},1}:\n [0.0]      \n [0.122174]\n [0.313727]\n [0.477212]\n [1.04692]  \n ⋮          \n [-0.636782]\n [-0.309111]\n [-0.154719]\n [-0.290622]\n [-0.290621]\n\njulia> plot(sol)\n","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: CAM1D Markov process time series)","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"One may also change the time step, dt, as well as the damping term, λ. Moreover, a seed may be specified in order to reproduce a time series over subsequent runs.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> sol = CAM1D(d, n, dt=1/12, λ=0.5, seed=42)\nretcode: Success\nInterpolation: 1st order linear\nt: 1001-element Array{Float64,1}:\n  0.0                \n  0.08333333333333333\n  0.16666666666666666\n  0.25               \n  0.3333333333333333\n  0.41666666666666663\n  0.49999999999999994\n  ⋮                  \n 82.83333333333299   \n 82.91666666666632   \n 82.99999999999964   \n 83.08333333333297   \n 83.1666666666663    \n 83.24999999999963   \n 83.25               \nu: 1001-element Array{Array{Float64,1},1}:\n [0.0]                 \n [-0.05538507236046933]\n [0.3759078513765658]  \n [0.31482891888314435]\n [0.2586995978771357]  \n [0.006314400534034209]\n [0.0819076515551696]  \n ⋮                     \n [0.6914165604620272]  \n [1.2048828846162347]  \n [0.6995798768325541]  \n [0.541755324779944]   \n [0.5939196502928672]  \n [0.6785936302181322]  \n [0.6785940957553153]  \n\njulia> plot(sol)\n","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: CAM1D Markov process time series with arguments)","category":"page"},{"location":"examples/#Hasselmann's-model-1","page":"Examples","title":"Hasselmann's model","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"It is also possible to create a Markov process time series using the Hasselmann (1976) climate model, described in:","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Hasselmann, K., 1976: Stochastic climate models Part I. Theory. Tellus, 28, 473–485, https://doi.org/10.1111/j.2153-3490.1976.tb00696.x.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"The Hasselmann Markov process – which is an AR(1) process – produces a time series that is normally distributed.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> n = 1000  # length of time series\n1000\n\njulia> sol = Hasselmann1D(n)\nretcode: Success\nInterpolation: 1st order linear\nt: 1001-element Array{Float64,1}:\n  0.0                 \n  0.041666666666666664\n  0.08333333333333333\n  0.125               \n  0.16666666666666666\n  0.20833333333333331\n  0.24999999999999997\n  ⋮                   \n 41.416666666666494   \n 41.45833333333316    \n 41.49999999999982    \n 41.54166666666649    \n 41.58333333333315    \n 41.624999999999815   \n 41.625               \nu: 1001-element Array{Array{Float64,1},1}:\n [0.0]                \n [0.02634147381570931]\n [0.06909671634959474]\n [0.311743726899033]  \n [0.5943736799301251]\n [0.4215357023451997]\n [0.2506514745054867]\n ⋮                    \n [-0.7733132505027923]\n [-0.8065075584943057]\n [-0.762325615339476]\n [-0.9280604718015151]\n [-1.2373636151238154]\n [-1.243929957136598]\n [-1.2439305063701647]\n\njulia> plot(sol)","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: Hasselmann1D Markov process time series)","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Again one may specify the time step, damping terms and the seed used for the random number generation.","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> sol = Hasselmann1D(n, dt=1/12, λ=0.5, seed=42)\nretcode: Success\nInterpolation: 1st order linear\nt: 1001-element Array{Float64,1}:\n  0.0                \n  0.08333333333333333\n  0.16666666666666666\n  0.25               \n  0.3333333333333333\n  0.41666666666666663\n  0.49999999999999994\n  ⋮                  \n 82.83333333333299   \n 82.91666666666632   \n 82.99999999999964   \n 83.08333333333297   \n 83.1666666666663    \n 83.24999999999963   \n 83.25               \nu: 1001-element Array{Array{Float64,1},1}:\n [0.0]                 \n [-0.13005535542166105]\n [0.04828669603525988]\n [0.6586911660509708]  \n [0.47664314043490363]\n [0.2527905702593072]  \n [0.5063085613562421]  \n ⋮                     \n [-0.470018895671266]  \n [-0.6640224731670596]\n [-0.6412697125734017]\n [-0.8226707386234937]\n [-0.9240097641730842]\n [-0.3550868399241891]\n [-0.3550865369728031]\n\njulia> plot(sol)\n","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"(Image: Hasselmann1D Markov process time series)","category":"page"},{"location":"examples/#Random-number-generation-1","page":"Examples","title":"Random number generation","text":"","category":"section"},{"location":"examples/#","page":"Examples","title":"Examples","text":"Produce a random number or an array of random numbers drawn from the SGS distribution d by invoking rand().","category":"page"},{"location":"examples/#","page":"Examples","title":"Examples","text":"julia> rand(d)\n-0.6568019443675626\n\njulia> rand(d, 10)\n10-element Array{Float64,1}:\n -0.6658576564565255\n -0.182584214456065  \n -0.3849235898350698\n -1.8058720924121987\n  0.08814146377788976\n  0.4741505815188713\n  0.20989773773209539\n  0.37743678587708196\n -0.19880370169821895\n  1.0225965534875916","category":"page"},{"location":"#SGSDist.jl-Documentation-1","page":"Home","title":"SGSDist.jl Documentation","text":"","category":"section"},{"location":"#Overview-1","page":"Home","title":"Overview","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"SGSDist.jl is a Julia package that implements the stochastically generated skewed (SGS) distribution, a probability density function (pdf) obtained from the correlated additive and multiplicative (CAM) noise stochastic climate model developed by Sura and Sardeshmukh (2008) and Sardeshmukh and Sura (2009).","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For more information on the form of the SGS distribution used in this package, see Sardeshmukh et al. (2015).","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"note: Note\nSGSDist.jl is not yet distributed through the Julia package manager, so one has to install the package locally.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"SGSDist.jl may be installed using the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"add SGSDist","category":"page"},{"location":"#Local-installation-1","page":"Home","title":"Local installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"If you want to install the package locally, first decide what directory will house the SGSDist.jl repository. There, clone the repository to your local machine:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"git clone git@github.com:brwst/SGSDist.jl.git","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Then, add the parent directory that contains the repository to the Julia LOAD_PATH by editing the startup.jl file located in ~/.julia/config:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"push!(LOAD_PATH, \"/path/to/dir/\")","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Reloading any Julia instances should make the SGSDist.jl module available for import.","category":"page"},{"location":"#References-1","page":"Home","title":"References","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"SGSDist.jl is based largely on the form of the SGS distribution and its associated Markov process outlined in Sardeshmukh et al. (2015):","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. Journal of Climate, 28, 9166–9187, https://doi.org/10.1175/JCLI-D-15-0020.1.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"For the theoretical development of the CAM noise stochastic climate model and the SGS distribution, see:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Sura, P., and P. D. Sardeshmukh, 2008: A Global View of Non-Gaussian SST Variability. Journal of Physical Oceanography, 38, 639–647, https://doi.org/10.1175/2007JPO3761.1.Sardeshmukh, P. D., and P. Sura, 2009: Reconciling Non-Gaussian Climate Statistics with Linear Dynamics. Journal of Climate, 22, 1193–1207, https://doi.org/10.1175/2008JCLI2358.1.","category":"page"}]
}