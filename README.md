# SGSDist.jl

`SGSDist.jl` is a Julia package that implements the stochastically generated skewed (SGS) distribution, a probability density function (pdf) obtained from the correlated additive and multiplicative (CAM) noise stochastic climate model developed by [Sura and Sardeshmukh (2008)](https://doi.org/10.1175/2007JPO3761.1) and [Sardeshmukh and Sura (2009)](https://doi.org/10.1175/2008JCLI2358.1).

For more information on the form of the SGS distribution used in this package, see [Sardeshmukh et al. (2015)](https://doi-org.proxy.lib.fsu.edu/10.1175/JCLI-D-15-0020.1).

## Installation

`SGSDist.jl` is not yet distributed through the Julia package manager, so one has to install the package locally.

First, clone the repository to your local machine:

```bash
git clone git@github.com:brwst/SGSDist.jl.git
```

Then, add the directory that contains the repository to the Julia `LOAD_PATH` by editing the `startup.jl` file located in `~/.julia/config`:

```
push!(LOAD_PATH, "/path/to/dir/")
```

You can then load the `SGSDist.jl` package with:

```julia
using SGSDist
```

## Documentation

More thorough documentation is found by building the documentation located in `doc` locally with the `Documenter.jl` package.

## Project Status

This package is being developed and tested against Julia 1.0 and above.
