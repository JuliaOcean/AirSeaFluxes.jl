var documenterSearchIndex = {"docs":
[{"location":"#AirSeaFluxes.jl-1","page":"Home","title":"AirSeaFluxes.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package is at a very early stage of development. Stay tuned ...","category":"page"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [AirSeaFluxes]","category":"page"},{"location":"#AirSeaFluxes.bulkformulae","page":"Home","title":"AirSeaFluxes.bulkformulae","text":"bulkformulae(atemp,aqh,speed,sst,hu=10,ht=2,hq=2,zref=10,atmrho=1.2)\n\nUnits: atemp  - mean air temperature (K)  at height ht (m) aqh    - mean air humidity (kg/kg) at height hq (m) speed  - mean wind speed (m/s)     at height hu (m) sst    - sea surface temperature (K)\n\nBulk formulae formulation:\n\nwind stress = (ust,vst) = rhoA * Cd * Ws * (del.u,del.v)\nSensib Heat flux = fsha = rhoA * Ch * Ws * del.T * CpAir\nLatent Heat flux = flha = rhoA * Ce * Ws * del.Q * Lvap\n                 = -Evap * Lvap\n\nwith Cd,Ch,Ce = drag coefficient, Stanton number and Dalton number respectively [no-units], function of height & stability; and\n\nWs = wind speed = sqrt(del.u^2 +del.v^2)\ndel.T = Tair - Tsurf ; del.Q = Qair - Qsurf\n\n\n\n\n\n","category":"function"},{"location":"#AirSeaFluxes.simpleflux-Tuple{Any,Any,Any}","page":"Home","title":"AirSeaFluxes.simpleflux","text":"simpleflux(Ca::Float,Co::Float,pisvel::Float)\n\nCompute flux entering the ocean per unit area\n\nmld=10 #mixed layer depth (m)\ntim=86400*30 #relaxation time scale (s)\npisvel=mld/tim #piston velocity (m/s)\nCo=0.0 #ocean value (e.g. concentration of some compound)\nCa=1.0 #atmospeheric value (e.g. equivalent compound concentration)\n\nsimpleflux(Ca,Co,pisvel)\n\n\n\n\n\n","category":"method"}]
}