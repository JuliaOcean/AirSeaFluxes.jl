var documenterSearchIndex = {"docs":
[{"location":"#AirSeaFluxes.jl-1","page":"Home","title":"AirSeaFluxes.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Computation and analysis of air-sea fluxes. ","category":"page"},{"location":"#","page":"Home","title":"Home","text":"warning: Warning\nThis package is in early developement stage when breaking changes can be expected.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [AirSeaFluxes]","category":"page"},{"location":"#AirSeaFluxes.AOGCM1D","page":"Home","title":"AirSeaFluxes.AOGCM1D","text":"AOGCM1D(ndays=60)\n\nAtmosphere-Ocean, coupled, one-dimensional column model.\n\noutputs,parameters=AOGCM1D(10)\np = dirname(pathof(AirSeaFluxes))\ninclude(joinpath(p,\"recipes_plots.jl\"))\np1,p2,p3=plot_final(outputs,parameters)\ndisplay(p3)\n\n\n\n\n\n","category":"function"},{"location":"#AirSeaFluxes.bulkformulae","page":"Home","title":"AirSeaFluxes.bulkformulae","text":"bulkformulae(atemp,aqh,speed,sst,hu=10,ht=2,hq=2,zref=10,atmrho=1.2)\n\nUnits: atemp  - mean air temperature (K)  at height ht (m) aqh    - mean air humidity (kg/kg) at height hq (m) speed  - mean wind speed (m/s)     at height hu (m) sst    - sea surface temperature (K)\n\nBulk formulae formulation:\n\nwind stress = (ust,vst) = rhoA * Cd * Ws * (del.u,del.v)\nSensib Heat flux = fsha = rhoA * Ch * Ws * del.T * CpAir\nLatent Heat flux = flha = rhoA * Ce * Ws * del.Q * Lvap\n                 = -Evap * Lvap\n\nwith Cd,Ch,Ce = drag coefficient, Stanton number and Dalton number respectively [no-units], function of height & stability; and\n\nWs = wind speed = sqrt(del.u^2 +del.v^2)\ndel.T = Tair - Tsurf ; del.Q = Qair - Qsurf\n\n\n\n\n\n","category":"function"},{"location":"#AirSeaFluxes.simpleflux-Tuple{Any, Any, Any}","page":"Home","title":"AirSeaFluxes.simpleflux","text":"simpleflux(Ca::Float,Co::Float,pisvel::Float)\n\nCompute flux entering the ocean per unit area\n\nmld=10 #mixed layer depth (m)\ntim=86400*30 #relaxation time scale (s)\npisvel=mld/tim #piston velocity (m/s)\nCo=0.0 #ocean value (e.g. concentration of some compound)\nCa=1.0 #atmospeheric value (e.g. equivalent compound concentration)\n\nsimpleflux(Ca,Co,pisvel)\n\n\n\n\n\n","category":"method"},{"location":"#AirSeaFluxes.angle_of_incidence-NTuple{4, Any}","page":"Home","title":"AirSeaFluxes.angle_of_incidence","text":"angle_of_incidence(lat,lon,jd,time)\n\nCalculate the cosine of the sun incident angle for a given location (lat,lon), julian day (jd) and time in hours and fraction of hours.\n\n\n\n\n\n","category":"method"},{"location":"#AirSeaFluxes.delta-Tuple{Any}","page":"Home","title":"AirSeaFluxes.delta","text":"delta(jd)\n\nCalculates the Sun's declination angle as a function of Julian day (jd)\n\n\n\n\n\n","category":"method"},{"location":"#AirSeaFluxes.heaviside-Tuple{Any}","page":"Home","title":"AirSeaFluxes.heaviside","text":"heaviside(x)\n\nHeaviside function : 0 when x<0 and 1 when x>=0\n\n\n\n\n\n","category":"method"}]
}
