function large(T,S,U,V,rho,Qnet,FWflux,ustar,ZF,lat)

#=
T=T(i,:)';
S=zeros(size(T));
U=(vlev:-1:1)/vlev;
V=zeros(size(U));
Qnet=Qnet(i,1);
FWflux=0;
ustar=ustar(i);
=#

#if ~exist('wmt','var')
#    load ltable.mat 'wmt' 'wst' 'usta' 'zehat';
#end

alpha=2e-4; #expansion coefficients for temperature
beta=7.4e-4; #expansion coefficients for salt
#tref=25; #reference temperature
#sref=25;
gravity_mks=9.81;
cp=4184; #J/K/kg
Rcrit=0.3;
karman=0.4;
rhow=1029; #kg/m^3
omega=7.2921e-5;

f=2*omega.*sin(lat/180*pi);

Z=ZF[2:end];
ZM=Z-diff(ZF)./2;
dz=-diff(ZF);
len=length(Z);

#rho=rhow.*(1 .-alpha.*(T.-tref).+beta.*(S.-sref));
#Br=gravity_mks.*(alpha.*tref+beta.*sref);
B=gravity_mks.*(alpha.*T+beta.*S);
Bf=gravity_mks.*(alpha.*Qnet     ./cp      ./rho[1]   .+beta.*FWflux./rho[1]);
# [m s-2]     [K-1]  [J s-1 m-2][J-1 kg K][kg-1 m^3]=
# [m s-2] [s-1 m] = [m^2 s-3]

L=ustar.^3 ./(karman.*Bf); #[m^3 s-3]*[s^3 m-2]=[m]
Cv=1.5; bt=0.2; Cs=98.96; epsln=0.1;
ws_star0=(abs.(Bf./Z)).^(1/3); #-sign(Bf).*
Vt2=max.(Cv.*sqrt.(bt)./(Rcrit.*karman.^2).*sqrt.(Cs.*epsln).*abs.(Z).*
    sqrt.(abs.([diff(B); 0]./dz)).*ws_star0,1e-8); #*ws

A=(0.1.*abs.(Z)).-abs.(Z').>0; #matrix with the indexes of surface layer
A[:,1]=max.(A[:,1],1);
N=sum(A,dims=2);

R=(B.-(sum(B'.*A,dims=2)./N)).*Z./
  ((sqrt.(U[1].^2+V[1].^2).-sqrt.(U.^2+V.^2)).^2 .+ Vt2.^2);
i=(R.>Rcrit).*[1:len;]; i=i[i.>0];
if isempty(i)
    h=abs(ZM[end]);
    i=length(R);
else
    h=abs(ZM[i[1]]);
    i=i[1];
end

if Bf>0
    he=0.7.*ustar./f;
    a=(abs.(Z).-he);
    tmp,ie=findmax(a.*(a.<0));
    a=(abs.(Z).-L);
    tmp,im=findmax(a.*(a.<0));
    i=min(ie[1],i);
    i=min(im[1],i);
end

i=max.(i,2);
imld=[1:i-1;];
ielse=[i:length(R);];

sigma=abs.(Z)./h;
xi=h/L;

i1=((sigma.>=0.1) .& (sigma.<1)); i2=sigma.<0.1;

wm=zeros(size(Z));
if !all(i1==0); wm[i1].=karman.*ustar./kpp_psim.(0.1.*xi); end
if !all(i2==0); wm[i2]=karman.*ustar./kpp_psim.(sigma[i2].*xi); end

ws=zeros(size(Z));
if !all(i1==0); ws[i1].=karman.*ustar./kpp_psis.(0.1.*xi); end
if !all(i2==0); ws[i2]=karman.*ustar./kpp_psis.(sigma[i2].*xi); end


ds=0.01;
sigmad=sigma.+ds;
i1=((sigmad.>=0.1) .& (sigmad.<1)); i2=sigmad.<0.1;

wmd=zeros(size(Z));
if !all(i1==0); wmd[i1].=karman.*ustar./kpp_psim.(0.1.*xi); end
if !all(i2==0); wmd[i2]=karman.*ustar./kpp_psim.(sigma[i2].*xi); end

wsd=zeros(size(Z));
if !all(i1==0); wsd[i1].=karman.*ustar./kpp_psis.(0.1.*xi); end
if !all(i2==0); wsd[i2]=karman.*ustar./kpp_psis.(sigma[i2].*xi); end

ni0=50e-4;
Ri0=0.7;
p1=3;

Rig=([-diff(B); 0]).*dz./[(diff(U).^2+diff(V).^2); 1];

nis=zeros(size(Rig));
nis[Rig.<0].=ni0;
nis[(Rig.>0) .& (Rig.<Ri0)]=ni0.*(1 .-(Rig[(Rig.>0) .& (Rig.<Ri0)]./Ri0).^2).^p1;
nis[Rig.>Ri0].=0;

nims=nis;
niss=nis;

nimw=1e-4;
nisw=1e-5;

nif=1e-3; Rrho0=1.9; p2=3;

Rrho=[alpha.*diff(T)./(beta.*diff(S)); 0];

nisd=zeros(size(Rrho));
nisd[(Rrho.>1) .& (Rrho.<Rrho0)]=nif.*(1 .-((Rrho[(Rrho.>1) .& (Rrho.<Rrho0)] .-1)./(Rrho0.-1)).^2).^p2;
nisd[(Rrho.<1) .| (Rrho.>Rrho0)].=0;

nitd=0.7.*nisd;

nim=nims.+nimw;
nis=niss.+nisw.+nisd;
nit=niss.+nisw.+nitd;

a1=1;
if length(imld)==length(R')
    nim0=nim[end];nis0=nis[end];nit0=nit[end];
    nim1=0;nis1=0;nit1=0;
else
    nim0=nim[imld[end]];nis0=nis[imld[end]];nit0=nit[imld[end]];
    nim1=nim[imld[end]+1];nis1=nis[imld[end]+1];nit1=nit[imld[end]+1];
end

#nim(i(1)+1)-nim(i(1))/(Z(i+1)-Z(i))

Gm1=nim0/(h.*wm[imld[end]]);
Gs1=nis0/(h.*ws[imld[end]]);
Gt1=nit0/(h.*ws[imld[end]]);
dsGm1=-(nim0-nim1)./dz[imld[end]]./wm[imld[end]]-nim0.*(wmd[imld[end]]-wm[imld[end]])./ds./(h.*wm[imld[end]].^2);
dsGs1=-(nis0-nis1)./dz[imld[end]]./ws[imld[end]]-nis0.*(wsd[imld[end]]-ws[imld[end]])./ds./(h.*ws[imld[end]].^2);
dsGt1=-(nit0-nit1)./dz[imld[end]]./ws[imld[end]]-nit0.*(wsd[imld[end]]-ws[imld[end]])./ds./(h.*ws[imld[end]].^2);

dsGm1=min(dsGm1,0);
dsGs1=min(dsGs1,0);
dsGt1=min(dsGt1,0);

a2m=-2+3 .*Gm1-dsGm1;
a2s=-2+3 .*Gs1-dsGs1;
a2t=-2+3 .*Gs1-dsGt1;

a3m=1-2 .*Gm1+dsGm1;
a3s=1-2 .*Gs1+dsGs1;
a3t=1-2 .*Gt1+dsGt1;

Gm=a1.*sigma+a2m.*sigma.^2+a3m.*sigma.^3;
Gs=a1.*sigma+a2s.*sigma.^2+a3s.*sigma.^3;
Gt=a1.*sigma+a2t.*sigma.^2+a3t.*sigma.^3;

Gm=(Gm[imld].*(Gm[imld].>0));
Gs=(Gs[imld].*(Gs[imld].>0));
Gt=(Gt[imld].*(Gt[imld].>0));

#=
Gm=sigma.*(1-sigma);
Gs=sigma.*(1-sigma);
Gt=sigma.*(1-sigma);
=#

kM=zeros(len);kT=zeros(len);kS=zeros(len);
kM[imld]=h.*wm[imld].*Gm[imld];
kT[imld]=h.*ws[imld].*Gt[imld];
kS[imld]=h.*ws[imld].*Gs[imld];

if !isempty(ielse)
    kM[ielse]=nim[ielse];
    kT[ielse]=nit[ielse];
    kS[ielse]=nis[ielse];
end

return kM,kT,kS

#dt=30;
#CFL=2.*dt./min(dz);

#kmax=1/CFL;

#kM=min(kM,kmax);
#kT=min(kT,kmax);
#kS=min(kS,kmax);

#[max(kT(imld)) max(ws(imld)) min(Gt(imld)) max(Gt(imld)) dsGt1 h 0.01.*Qnet]

end
