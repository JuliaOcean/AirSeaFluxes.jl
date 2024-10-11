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

    atemp=300.
    aqh=0.001
    speed=1.
    sst=10.
    all=bulkformulae(atemp,aqh,speed,sst)
    @test isapprox(all.hl,-3.0606099804357885; rtol=1e-2)
    @test isapprox(all.hs,2.0282408526727473; rtol=1e-2)
    @test isapprox(all.evap,1.2244888899523058e-9; rtol=1e-2)
    @test isapprox(all.tau,0.00915548218468587; rtol=1e-2)

    outputs,parameters=AOGCM1D(10)
    @test isa(outputs,Dict)

    pco2=calc_pco2(10,35,2.1,2.3)
    @test isapprox(pco2,0.5393241959361947; rtol=1e-2)
end
