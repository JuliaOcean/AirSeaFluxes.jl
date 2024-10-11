using Documenter, AirSeaFluxes

makedocs(;
    modules=[AirSeaFluxes],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaOcean/AirSeaFluxes.jl/blob/{commit}{path}#L{line}",
    sitename="AirSeaFluxes.jl",
    authors="JuliaOcean <gforget@mit.edu>",
    warnonly = [:cross_references,:missing_docs],
)

deploydocs(;
    repo="github.com/JuliaOcean/AirSeaFluxes.jl",
)
