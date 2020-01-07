# SGSDist.jl Documentation

## Overview

`SGSDist.jl` is a Julia package that implements the stochastically generated skewed (SGS) distribution, a probability density function (pdf) obtained from the correlated additive and multiplicative (CAM) noise stochastic climate model developed by [Sura and Sardeshmukh (2008)](https://doi.org/10.1175/2007JPO3761.1) and [Sardeshmukh and Sura (2009)](https://doi.org/10.1175/2008JCLI2358.1).

For more information on the form of the SGS distribution used in this package, see [Sardeshmukh et al. (2015)](https://doi-org.proxy.lib.fsu.edu/10.1175/JCLI-D-15-0020.1).

## Installation

!!! note

    `SGSDist.jl` is not yet distributed through the Julia package manager, so one has to install the package locally.

`SGSDist.jl` may be installed using the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```julia
add SGSDist
```

### Local installation

If you want to install the package locally, first decide what directory will house the `SGSDist.jl` repository. There, clone the repository to your local machine:

```
git clone git@github.com:brwst/SGSDist.jl.git
```

Then, add the parent directory that contains the repository to the Julia `LOAD_PATH` by editing the `startup.jl` file located in `~/.julia/config`:

```julia
push!(LOAD_PATH, "/path/to/dir/")
```

Reloading any Julia instances should make the `SGSDist.jl` module available for import.

## References

`SGSDist.jl` is based largely on the form of the SGS distribution and its associated Markov process outlined in [Sardeshmukh et al. (2015)](https://doi-org.proxy.lib.fsu.edu/10.1175/JCLI-D-15-0020.1):

>Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015: Need for Caution in Interpreting Extreme Weather Statistics. *Journal of Climate*, **28**, 9166–9187, [https://doi.org/10.1175/JCLI-D-15-0020.1](https://doi.org/10.1175/JCLI-D-15-0020.1).

For the theoretical development of the CAM noise stochastic climate model and the SGS distribution, see:

>Sura, P., and P. D. Sardeshmukh, 2008: A Global View of Non-Gaussian SST Variability. *Journal of Physical Oceanography*, **38**, 639–647, [https://doi.org/10.1175/2007JPO3761.1](https://doi.org/10.1175/2007JPO3761.1).
>
>Sardeshmukh, P. D., and P. Sura, 2009: Reconciling Non-Gaussian Climate Statistics with Linear Dynamics. *Journal of Climate*, **22**, 1193–1207, [https://doi.org/10.1175/2008JCLI2358.1](https://doi.org/10.1175/2008JCLI2358.1).
