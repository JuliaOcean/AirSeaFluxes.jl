function holtslag(Tv,q,U,V,LH,SH,ustar,ZF)
# a vertical diffusion scheme based on Holtslag and Boville (1993)
# Includes enhenced diffusion below the PBL

#TF-turbulent fluxes, positive upward
#Tv-virtual potential temperature

Rcrit=0.5;
gravity_mks=9.81;
karman=0.4;
a=7.8;
eps=0.1;
rhoa=1.2;
cpa=1005; #J/K/kg
Av=2500000;
lamdac=300;
humid_fac=0.606;

Z=ZF[2:end];
DZ=diff(ZF,dims=1);
ZM=Z-DZ./2;
SZM=length(ZM);

ZD=Z[1]:0.01*(Z[1]):Z[end];

#TF=SH.*(1 .+humid_fac.*q)./rhoa./cpa.+humid_fac.*(Tv[1]./(1 .+humid_fac*q)).*LH./rhoa./Av;
TF=SH./rhoa./cpa.+humid_fac.*Tv[1].*LH./rhoa./Av;

#TF=(SH+LH)./rhoa./cpa;
ustar=abs(ustar)
L=Tv[1].*ustar.^3 ./(karman.*gravity_mks.*TF);
HD=eps.*ZD./L;
H=eps.*ZF[2:end]./L;

wsd=ustar./holtslag_psis.(HD);
Ts=Tv[1].+max.(min.(a.*TF./wsd,3),0);

spl = Spline1D(ZM, U); UD=spl(ZD);
spl = Spline1D(ZM, V);  VD=spl(ZD);
spl = Spline1D(ZM, Tv); TD=spl(ZD);

ha=abs.(Rcrit.*Ts.*(UD.^2+VD.^2)./(gravity_mks.*(TD-Ts)))-ZD;
ha[ha.>0].=Inf;
if all(isinf.(ha))
    if TF>0
        id=SZM;
    else
        id=1;
    end
else
    ~,id=findmin(abs.(ha));
end
id=max(id,1);
h=abs(ZD[id[1]]);
n=(h.-ZF).*((h.-ZF).>0);
n[n.==0].=Inf;
~,i=findmin(n);
frac=(h-ZF[i])./DZ[i];

nsl=(0.1 .*h.-ZF).*((0.1 .*h.-ZF).>0);
nsl[nsl.==0].=Inf;
~,isl=findmin(nsl);

ipbl=[1:i;];
if i==length(ZM)
    npbl=[];
else
    npbl=i+1:length(ZM);
end

kM=zeros(SZM);
kS=zeros(SZM);
wm=zeros(SZM);
ws=((gravity_mks/Tv[1]).*abs(TF).*h).^(1/3);

wm[1:isl]=ustar./holtslag_psim.(H[1:isl]);
wm[isl+1:i].=(ustar.^3 .+0.6.*ws.^3).^(1/3);

kM[ipbl]=karman.*(ZF[ipbl]+DZ[ipbl]/2).*
     wm[i].*(1 .-(ZF[ipbl]+DZ[ipbl]/2)./h).^2;
Pr=(holtslag_psis.(H[isl])./
    holtslag_psim.(H[isl])+
    a.*karman.*0.1.*ws./wm[isl]);
kS[ipbl]=kM[ipbl]./Pr;

S=abs.(diff([0; U])./diff(ZF));
GRi=(gravity_mks./Tv).*diff([Tv[1]; Tv])./diff(ZF)./S.^2;

fm=zeros(SZM);
fm[GRi.>=0]=1 ./(1 .+10 .*GRi[GRi.>=0].*(1 .+8 .*GRi[GRi.>=0]));
fm[GRi.<0]=(1 .-18 .*GRi[GRi.<0]).^(0.5);

lc=1 ./((1 ./(karman.*diff(ZF))).+(1 ./(karman.*lamdac)));

kM[npbl]=lc[npbl].^2 .*S[npbl].*fm[npbl];
kS[npbl]=kM[npbl];

kM[i]=kM[i].*frac+lc[i].^2 .*S[i].*fm[i].*(1-frac);
kS[i]=kS[i].*frac+lc[i].^2 .*S[i].*fm[i].*(1-frac);

return kM,kS
#[i max(kM) max(kS) TF max(Ts) U(1)]

end
