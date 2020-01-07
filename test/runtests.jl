using Pkg
using Random
using Test

tests = []

if length(ARGS) > 0
    tests = ARGS
end

Random.seed!(42)    # set a prng seed
tolerance = 1e-8    # set error tolerance

for t in tests
    include(t)      # include the different test suites
end
