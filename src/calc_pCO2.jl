# calculate pCO2 for given Temperature(T), Salinity(S), DIC, Alkalinity, etc..
# Efficient solver following Follows et al (2005)
# pco2eq = atmospheric reference pCO2 level (atmospheres)
#                for which to find equilibrium dic, csat
#       csat = equilibrium total inorganic carbon (mol/m³)
#             where 1 T = 1 metric ton = 1000 kg
#       ta  = total alkalinity (eq/m³)
#       pt  = inorganic phosphate (mol/m³)
#       sit = inorganic silicate (mol/m³)
#       T   = temperature (degrees C)
#       S   = salinity (PSU)
#		    hg  = first guess of [H+] (10e-8 is a good cold start value)

# set pt = 0.0, sit = 0.0, if the model doesn't have these two values
function calc_pco2(T,S,dic,ta,pt=0.0,sit=0.0)
    # set first guess for [H⁺]
    hg = 1.0e-8
    # convert degrees C to degrees K
    TK = T + 273.0

    # estimated concentration of borate(bt) based on salinity
    scl=S/1.80655;
    bt=0.000232 * scl/10.811

    # Coefficient algorithms as used in OCMIP2 protocols K1, K2 Millero (1995) using Mehrbach data
    k1 = 10^(-1.0*(3670.7/TK - 62.008 + 9.7944*log(TK) - 1.18e-2*S + 1.16e-4*S*S))
    k2 = 10^(-1.0*(1394.7/TK + 4.777 - 1.84e-2*S + 1.18e-4*S*S))

    # K1p, K2p, K3p, DOE (1994)
    k1p = exp(-4576.752/TK + 115.525 - 18.453*log(TK) + (-106.736/TK+ 0.69171)*sqrt(S)+ (-0.65643/TK- 0.01844)*S)
    k2p = exp(-8814.715/TK + 172.0883- 27.927*log(TK) + (-160.34/TK + 1.3566)*sqrt(S) + (0.37335/TK - 0.05778)*S)
    k3p = exp(-3070.75/TK + 18.141 + (17.27039/TK + 2.81197)*sqrt(S) + (-44.99486/TK - 0.09984)*S)

    # Kb, Millero (1995) using data from Dickson
    kb = exp((-8966.90 - 2890.53*sqrt(S) - 77.942*S + 1.728*S^1.5 - 0.0996*S*S)/TK + (148.0248 + 137.1942*sqrt(S) + 1.62142*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S) * log(TK) + 0.053105*TK*sqrt(S))

    # Kw, Millero (1995)
    kw = exp(-13847.26/TK + 148.9652 - 23.6521*log(TK) + (118.67/TK - 5.977 + 1.0495*log(TK))*sqrt(S) - 0.01615*S)

    # fugacity, Weiss and Price, Marine Chem, 8, 347 (1990)
    TK1 = TK/100.0
    ff = exp(-162.8301 + 218.2968/TK1 + 90.9241*log(TK1) - 1.47696*(TK1*TK1) + S*(0.025695 - 0.025225*TK1 + 0.0049867*(TK1*TK1)))

    # Ksi, Millero (1995) if needed
    I = (19.924*S)/(1000 - 1.005*S)
    ksi = exp(-8904.2/TK + 117.385 - 19.334*log(TK) + (-458.79/TK + 3.5913)*sqrt(I) + (188.74/TK - 1.5998)*I + (-12.1652/TK + 0.07871)*I*I + log(1.0 - 0.001005*S))

    # First guess of [H⁺]: from last timestep *OR* fixed for cold start here iterate for accurate solution
    for i in 1:1000
        # estimate contributions to total alk from borate, silicate, phosphate
        bohg = (bt*kb)/(hg+kb)
        siooh3g = (sit*ksi)/(ksi + hg)
        denom = (hg*hg*hg) + (k1p*hg*hg) + (k1p*k2p*hg) + (k1p*k2p*k3p)
        h3po4g = (pt*hg*hg*hg)/denom
        h2po4g = (pt*k1p*hg*hg)/denom
        hpo4g = (pt*k1p*k2p*hg)/denom
        po4g = (pt*k1p*k2p*k3p)/denom
        # estimate carbonate alkalinity
        cag = ta - bohg - (kw/hg) + hg - hpo4g - 2*po4g + h3po4g - siooh3g
        # improced estimate of hydrogen ion concentration
        gamm = dic/cag
        dummy = (1 - gamm)*(1 - gamm)*k1*k1 - 4*k1*k2*(1 - 2*gamm)
        Hnew = 0.5*((gamm - 1)*k1 + sqrt(dummy))
        hg = Hnew
    end
    co2s = dic/(1 + (k1/hg) + ((k1*k2)/(hg*hg)))
    pco2 = co2s/ff
    return pco2
end
