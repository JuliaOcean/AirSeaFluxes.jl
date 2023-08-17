"""
    simpleflux(Ca::Float,Co::Float,pisvel::Float)

Compute flux entering the ocean per unit area

```
mld=10 #mixed layer depth (m)
tim=86400*30 #relaxation time scale (s)
pisvel=mld/tim #piston velocity (m/s)
Co=0.0 #ocean value (e.g. concentration of some compound)
Ca=1.0 #atmospeheric value (e.g. equivalent compound concentration)

simpleflux(Ca,Co,pisvel)
```
"""
@inline simpleflux(Ca,Co,pisvel) = pisvel*(Ca-Co)
