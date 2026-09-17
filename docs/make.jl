using Documenter
using TurbulentAnalysis

makedocs(
    sitename = "TurbulentAnalysis",
    format = Documenter.HTML(),
    modules = [TurbulentAnalysis],
    remotes = nothing,
    pages = ["Library" => "library.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
