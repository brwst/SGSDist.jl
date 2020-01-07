"""

    air()

Produce an array of sample climate data for testing. This sample time
series is taken from the [NCEP/NCAR Reanalysis 1 dataset](https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2),
representing the 925 hPa air temperature standardized anomalies of all 
DJF winters from 1948-2018 at 50.0°N, 235.0°E.

See also:

>Kalnay, E., and Coauthors, 1996: The NCEP/NCAR 40-year reanalysis project. *Bulletin of the American Meteorological Society*, **77**, 437–471, [https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2](https://doi.org/10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2).


# Examples
```jldoctest
julia> air()
6408-element Array{Float64,1}:
  0.67547214
  0.8308915
  0.19012427
 -0.20711477
 -0.10181273
  0.7421673
  0.22026125
  0.050515447
  0.043254483
  0.044615462
  ⋮
  0.7591071
  0.54696834
  0.3421846
  0.044188347
  0.2433472
  0.4460243
  1.4941149
  0.49388295
 -0.052520715

```

"""
function air()
    return vec(readdlm(joinpath(dirname(pathof(SGSDist)),
                                "../test/data/R1_air_stdanom_t1948-2018-djf_x235p0_y50p0_z925.csv"),' ','\n'))
end
