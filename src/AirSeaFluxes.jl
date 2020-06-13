module AirSeaFluxes

using Dierckx
using LinearAlgebra

export simpleflux, bulkformulae, AOGCM1D

include("simpleflux.jl")
include("bulkformulae.jl")
include("AOGCM1D.jl")

end # module
