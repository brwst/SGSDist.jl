# SGSDist.jl

| **Documentation**                 | **Build Status**                                                                                |
|:---------------------------------:|:------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][travis-img]][travis-url]  |

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

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**
- [**DEV**][docs-dev-url] &mdash; *documentation of the in-development version.*


[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://brwst.github.io/SGSDist.jl/dev

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://brwst.github.io/SGSDist.jl/stable

[travis-img]: https://travis-ci.org/brwst/SGSDist.jl.svg?branch=master
[travis-url]: https://travis-ci.org/brwst/SGSDist.jl
