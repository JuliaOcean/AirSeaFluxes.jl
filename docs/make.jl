using Documenter, AirSeaFluxes

makedocs(;
    modules=[AirSeaFluxes],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gaelforget/AirSeaFluxes.jl/blob/{commit}{path}#L{line}",
    sitename="AirSeaFluxes.jl",
    authors="gaelforget <gforget@mit.edu>",
    assets=String[],
)

deploydocs(;
    repo="github.com/gaelforget/AirSeaFluxes.jl",
)
