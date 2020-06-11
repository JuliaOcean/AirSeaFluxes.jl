function bulk(atemp,aqh,speed,sst,hu=10,ht=2,hq=2,zref=10,atmrho=1.2)
# formulae in short:
#   wind stress = (ust,vst) = rhoA * Cd * Ws * (del.u,del.v)
#   Sensib Heat flux = fsha = rhoA * Ch * Ws * del.T * CpAir
#   Latent Heat flux = flha = rhoA * Ce * Ws * del.Q * Lvap
#                    = -Evap * Lvap
# with Ws = wind speed = sqrt(del.u^2 +del.v^2) ;
#      del.T = Tair - Tsurf ; del.Q = Qair - Qsurf ;
#      Cd,Ch,Ce = drag coefficient, Stanton number and Dalton number
#            respectively [no-units], function of height & stability

umin=0.1;
karman=0.4;
gravity_mks=9.81;
cen2kel=273.150;
humid_fac=0.606;
cvapor_fac=640380;
#cvapor_fac_ice=11637800;
cvapor_exp=5107.4;
#cvapor_exp_ice=5897.8;
saltsat=0.980;
#atmrho=1.2;
gamma_blk=0.01;
cdrag_1=0.0027000;
cdrag_2=0.0001420;
cdrag_3=0.0000764;
cstanton_1 = 0.0327;
cstanton_2 = 0.0180;
cdalton = 0.0346;
niter_bulk=2;
psim_fac=5;
atmcp=1005;
flamb=2500000;
rhoConstFresh=999.8;
stefanBoltzmann = 5.670e-8;
ocean_emissivity=5.50e-8 / 5.670e-8;
albedo=0.1;

#-- Set surface parameters :
zwln = log(hu./zref);
ztln = log(ht./zref);
czol = hu*karman*gravity_mks;

wspeed=max(speed,umin);

#-   Surface Temp.
Tsf = sst + cen2kel;

#--- Compute turbulent surface fluxes
#-   Pot. Temp and saturated specific humidity

t0 = atemp.*(1 + humid_fac*aqh);
tmpbulk = cvapor_fac*exp(-cvapor_exp./Tsf);
ssq = saltsat*tmpbulk/atmrho;
deltap = atemp + gamma_blk*ht - Tsf;
delq   = aqh - ssq;

#--  initial guess for exchange coefficients:
#    take U_N = del.U ; stability from del.Theta ;
stable = 0.5 + 0.5*sign(deltap);

tmpbulk = cdrag_1./wspeed + cdrag_2 + cdrag_3.*wspeed; # CD

rdn = sqrt(tmpbulk);
rhn = (1-stable)*cstanton_1 + stable*cstanton_2;
ren = cdalton;
#--  calculate turbulent scales
ustar=rdn.*wspeed;
tstar=rhn.*deltap;
qstar=ren*delq;

rd=0;rh=0;re=0;huol=0;psimh=0;psixh=0;

#--- iterate with psi-functions to find transfer coefficients
for iter=1:niter_bulk

    huol = ( tstar./t0 + qstar./(1/humid_fac+aqh) ).*czol./(ustar.*ustar);

    #cph The following is different in Large&Pond1981 code:
    #cph huol = max(huol,zolmin) with zolmin = -100
    tmpbulk = min(abs(huol),10);
    huol   = tmpbulk.*sign(huol);
    htol   = huol.*ht./hu;
    hqol   = huol.*hq./hu;
    stable = 0.5 + 0.5*sign(huol);

    #                 Evaluate all stability functions assuming hq = ht.
    # The following is different in Large&Pond1981 code:
    # xsq = max(sqrt(abs(1 - 16.*huol)),1)
    xsq    = sqrt( abs(1 - huol*16) );
    x      = sqrt(xsq);
    psimh = -psim_fac*huol.*stable + (1-stable).*
        ( log( (1 + 2*x + xsq).*(1+xsq)*.125 ) - 2*atan(x) + 0.5*pi );

    # The following is different in Large&Pond1981 code:
    # xsq = max(sqrt(abs(1 - 16.*htol)),1)
    xsq   = sqrt( abs(1 - htol*16) );
    psixh = -psim_fac*htol.*stable + (1-stable).*( 2*log(0.5*(1+xsq)) );

    #-   Shift wind speed using old coefficient
    usn = wspeed./(1 + rdn/karman.*(zwln-psimh) );
    usm = max(usn, umin);

    #-   Update the 10m, neutral stability transfer coefficients
    tmpbulk = cdrag_1./usm + cdrag_2 + cdrag_3.*usm;

    rdn = sqrt(tmpbulk);
    rhn = (1-stable)*cstanton_1 + stable*cstanton_2;
    ren = cdalton;

    #-   Shift all coefficients to the measurement height and stability.
    rd = rdn./(1 + rdn.*(zwln-psimh)/karman);
    rh = rhn./(1 + rhn.*(ztln-psixh)/karman);
    re = ren./(1 + ren.*(ztln-psixh)/karman);

    #--  Update ustar, tstar, qstar using updated, shifted coefficients.
    ustar = rd.*wspeed;
    qstar = re.*delq;
    tstar = rh.*deltap;

    # end of iteration loop
end

#-   Coeff:
tau   = atmrho*rd.*wspeed;

#-   Turbulent Fluxes
hs = atmcp*tau.*tstar;
hl = flamb*tau.*qstar;
#   change sign and convert from kg/m^2/s to m/s via rhoConstFresh
evap = -tau.*qstar/rhoConstFresh;
ch=rd*rh
ce=rd*re

return hl,hs,evap,ch,ce,tau,ssq,huol,rd,re,rh,ustar,qstar,tstar,psimh,psixh

end
