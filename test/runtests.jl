using AirSeaFluxes
using Test

@testset "AirSeaFluxes.jl" begin
    mld=10 #mixed layer depth (m)
    tim=86400*30 #relaxation time scale (s)
    pisvel=mld/tim #piston velocity (m/s)
    Co=0.0 #ocean value (e.g. concentration of some compound)
    Ca=1.0 #atmospeheric value (e.g. equivalent compound concentration)
    flx=simpleflux(Ca,Co,pisvel)
    @test isapprox(flx,3.858024691358025e-6; atol=1e-2)
end
