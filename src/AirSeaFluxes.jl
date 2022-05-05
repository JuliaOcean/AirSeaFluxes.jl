module AirSeaFluxes

using Dierckx
using LinearAlgebra

export simpleflux, bulkformulae, AOGCM1D

include("simpleflux.jl")
include("bulkformulae.jl")
include("holtslag.jl")
include("kpp.jl")
include("AOGCM1D_helpers.jl")
include("AOGCM1D.jl")
include("calc_pCO2.jl")

end # module
