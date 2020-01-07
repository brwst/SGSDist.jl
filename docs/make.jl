using SGSDist
using Documenter

DocMeta.setdocmeta!(SGSDist, :DocTestSetup, :(using SGSDist); recursive=true)

makedocs(
    sitename = "SGSDist.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Robert West",
    modules = [SGSDist],
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Library" => "library.md",
        "Index" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/brwst/SGSDist.jl.git",
)
