using Dierckx
using Plots
using RollingFunctions
using LinearAlgebra

#using PyPlot
include("heaviside.jl")
include("delta.jl")
include("angle_of_incidence.jl")
include("kpp_psim.jl")
include("kpp_psis.jl")
include("bulk.jl")
include("holtslag.jl")
include("holtslag_psim.jl")
include("holtslag_psis.jl")
include("large.jl")

cpl=1;
moist=1;

dt=60; #seconds
ndays=60; # duration of simulation in days

n=Int(round(ndays*24*60*60/dt)); # number of time steps

H=50; # atmospheric layer thickness
dz=1; # oceanic layer thickness
vlevo=50; # number of oceanic vertical levels
vleva=50; # number of atmospheric vertical levels
rhow0=1029; # kg/m^3
cp=4184; #J/K/kg
cpa=1005; #J/K/kg
sb = 5.670e-8;
d2k=273.15;
Av=2500000;
gamma_blk=0.01;
cvapor_fac=640380;
cvapor_exp=5107.4;
cvapor_exp_ice=5897.8;
saltsat=0.980;
p0=101300; #reference pressure
gravity_mks=9.81;
humid_fac=0.606;
tref=25; #reference temperature
sref=25;
SO=sref; #salinity not implemented
alpha=2e-4; #expansion coefficients for temperature
beta=7.4e-4; #expansion coefficients for salt
emissivity=0.97; # ocean emissivity
kappa=1e-5; # m^2/kg air absorption coefficient
Rgas=287.05;
karman=0.4;

ZA=H*[1:vleva;];
ZO=zeros(vlevo+1);
ZO[2:vlevo+1]=-[1:vlevo;]*dz;
kT0=ones(vlevo)*1e-4; # SST vertical mixing
kT=kT0; # SST vertical mixing

Aimp=0.5; # Atmospheric implicitness
Oimp=0.5; # Oceanic implicitness

R=0.62; g1=0.6; g2=20;

f=R*exp.(ZO/g1)+(1-R)*exp.(ZO/g2);
f[ZO.<-200].=0;
f=-diff(f,dims=1);

albedo=0.06;
lat=30;
lon=0;
t=[1:n;]*dt/3600;
jstart=70;# julian day, borial summer
jadd=round.((t.-t[1])/24);
jd=jstart.+jadd;

AOI=angle_of_incidence.(lat,lon,jd,t); AOI[AOI.<0].=0;
SW=(1362*(1-albedo).*f'.*AOI);
#mean(SW(:,1))
#SW=(repmat((sin(2*pi/n*(1:n)*ndays+pi/2)+1)/2*1362*(1-albedo),vlev,1).*
#repmat(f',1,n))'.*cos(lat/180*pi);

TO=zeros(n+1,vlevo);
UO=zeros(n+1,vlevo);
SST=zeros(n+1);

UA=zeros(n+1,vleva);
qA=zeros(n+1,vleva);
THETA=zeros(n+1,vleva);
rhoa=zeros(n+1,vleva);
TA=zeros(n+1,vleva)

UOst=zeros(vlevo); TOst=zeros(vlevo);
UAst=zeros(vleva); qAst=zeros(vleva); THETAst=zeros(vleva);

############################################################################
# Initial conditions
############################################################################
#U0=ones(n,1)*5;
#THETA0=(sin(2*pi/n*(1:n)*ndays/5-pi/2)+21);
TA[1,:]=25 .-0.007.*H.*[1:vleva;];
P=zeros(vleva);
P[1]=101300;
P[2:vleva]=P[1].*exp.(-cumsum(H./(Rgas.*(TA[1,2:end].+d2k)./gravity_mks)));
rhoa[1,:]=P./(Rgas.*(TA[1,:].+d2k));
THETA[1,:]=(TA[1,:].+d2k).*(P[1]./P).^(2/7).-d2k;
THETA[1,1:20].=THETA[1,1];
THETA0=ones(n,1)*THETA[1,end].+1;#+gamma_blk*H+0.01;#+0.51;#+SW*kappa/cpa*dt;
T0=ones(n+1,vlevo).*26;#repmat((16+(vlev-1:-1:0)/vlev*8),n+1,1);
TO[:,:].=T0;
TO[:,1]=T0[:,1]; TO[1,:]=T0[1,1].+15 .*([vlevo:-1:1;]/vlevo.-1);
SST.=TO[:,1];
ssq = saltsat*cvapor_fac*exp(-cvapor_exp./(TA[1]+d2k))./rhoa[1,1];
RH=0.8;
q0=0.5*ssq;
#Tv=(THETA(i-1,:)+d2k).*(1+humid_fac*q(i-1));
qA[1,:]=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA[1,:].+d2k))./rhoa[1,:].*RH;
UA[1,:]=(0.2/karman).*(log.((ZA.-H/2)./1e-4));
#U(1,:)=U(1,end);
U0=(0.2/karman).*(log((ZA[end]+H/2)/1e-4)).*ones(1,n+1);
atemp=THETA[1,1]+d2k; aqh=qA[1,1];
speed=abs(UA[1,1]); sst=SST[1];
tmp0=bulk(atemp,aqh,speed,sst,H/2,H/2,H/2,H/2,rhoa[1,1]);
huol0=tmp0[8];
TAup=270;
#U0=(sin(2*pi/n*(1:n)*ndays/5+2*pi/n*3600/dt*sf)*8+U(1,end));
#UO(:)=0.;
#U0(3600*24*5/dt:end)=U(1,end); one pertubation
############################################################################

LH=zeros(1,n); SH=zeros(1,n); TAU=zeros(1,n); TAUO=zeros(1,n);
LW=zeros(1,n); E=zeros(1,n); ssq=zeros(1,n);
xi=zeros(1,n); ce=zeros(1,n); ch=zeros(1,n);
ustar=zeros(1,n); qstar=zeros(1,n); tstar=zeros(1,n);
psimh=zeros(1,n); psixh=zeros(1,n); Ri=zeros(n,vleva);
rd=zeros(1,n); rh=zeros(1,n); re=zeros(1,n);
Qnet=zeros(1,n);
KAm=zeros(n,vleva);KAt=zeros(n,vleva);
KOm=zeros(n,vlevo);KOt=zeros(n,vlevo);

for i=2:n+1
    println(i)
    if rem(i,1440*60/dt)==0
        println(string("Day  ",round(i/(1440*60/dt),digits=5)))
    end

    ########################################################################
    # surface
    ########################################################################
    hl,hs,_,ce[i-1],ch[i-1],_,ssq[i-1],xi[i-1],rd[i-1],
    rh[i-1],re[i-1],ustar[i-1],qstar[i-1],tstar[i-1],psimh[i-1],psixh[i-1]=
    bulk(TA[i-1,1]+d2k,qA[i-1,1],abs(UA[i-1,1]-UO[i-1,1]),SST[i-1,1],
        H/2,H/2,H/2,10,rhoa[i-1,1]);
    LH[i-1]=-hl;
    SH[i-1]=-hs;
    E[i-1]=-hl/Av;
    ϵa=0.725+0.17*log10(sum(qA[1,:]/1000*1.2*H)*100)
    LW[i-1]=emissivity*sb*((SST[i-1]+d2k).^4 .- ϵa.*(TAup).^4);
    TAU[i-1]=rhoa[i-1,1].*ustar[i-1].^2 .*sign(ustar[i-1]);
    Qnet[i-1]=SW[i-1,1]-LH[i-1]-SH[i-1]-LW[i-1];
    ########################################################################
    # Atmospheric vertical diffusion
    ########################################################################

    Tv=(THETA[i-1,:].+d2k).*(1 .+humid_fac*qA[i-1,:]);
    KAm[i-1,:], KAt[i-1,:]= holtslag(Tv,qA[i-1,1],UA[i-1,:],
                 UA[i-1,:].*0,LH[i-1],SH[i-1],ustar[i-1],[0; ZA]);
    kam=KAm[i-1,:];
    kat=KAt[i-1,:];

    ########################################################################
    # Oceanic vertical diffusion
    ########################################################################

    rhow=rhow0.*(1 .-alpha.*(TO[i-1,:].-tref).+beta.*(SO.-sref));
    B=gravity_mks.*(rhow0.-rhow)./rhow0;
    Uo=(vlevo:-1:1)/vlevo/10; Vo=zeros(size(Uo)); FWflux=zeros(size(Qnet));
    KOm[i-1,:],KOt[i-1,:],_=
        large(TO[i-1,:],TO[i-1,:].*0,UO[i,:],Vo,rhow,Qnet[i-1]
             ,0,ustar[i-1],ZO,lat);
    kot=KOt[i-1,:];
    kom=KOm[i-1,:];

    ########################################################################
    # Atmospheric explicit time stepping
    ########################################################################

    TAUS=sign(UA[i-1,1]-UO[i-1,1]);

    UAst[1]=UA[i-1,1]-TAU[i-1]*TAUS/(H*rhoa[i-1,1])*dt-(1-Aimp)*(
        kam[1].*(UA[i-1,1]-UA[i-1,2])/H^2*dt);
    UAst[end]=UA[i-1,end]-(1-Aimp)*(
        kam[end]*(UA[i-1,end]-U0[i-1])/H^2*dt);
    UAst[2:end-1]=UA[i-1,2:end-1]+(1-Aimp)*(
        kam[1:end-2].*(UA[i-1,1:end-2]-UA[i-1,2:end-1])/H^2*dt-
        kam[2:end-1].*(UA[i-1,2:end-1]-UA[i-1,3:end  ])/H^2*dt);

    qAst[1]=qA[i-1,1]-E[i-1]/(H*rhow[1])*dt-
        (1-Aimp)*(kat[1]*(qA[i-1,1]-qA[i-1,2])/H^2*dt);
    qAst[end]=qA[i-1,end]-(1-Aimp)*(kat[end]*(qA[i-1,end]-q0)/H^2*dt);
    qAst[2:end-1]=qA[i-1,2:end-1]+(1-Aimp)*(
        kat[1:end-2].*(qA[i-1,1:end-2]-qA[i-1,2:end-1])/H^2*dt-
        kat[2:end-1].*(qA[i-1,2:end-1]-qA[i-1,3:end  ])/H^2*dt);

    THETAst[1]=THETA[i-1,1]+(SH[i-1])/cpa/(H*rhoa[i-1,1])*dt-
         (1-Aimp)*(kat[1]*(THETA[i-1,1]-THETA[i-1,2])/H^2*dt);
    THETAst[end]=THETA[i-1,end]-
         (1-Aimp)*(kat[end]*(THETA[i-1,end]-THETA0[i-1])/H^2*dt);
    THETAst[2:end-1]=THETA[i-1,2:end-1]+(1-Aimp)*(
         kat[1:end-2].*(THETA[i-1,1:end-2]-THETA[i-1,2:end-1])/H^2*dt-
         kat[2:end-1].*(THETA[i-1,2:end-1]-THETA[i-1,3:end  ])/H^2*dt);

    ########################################################################
    # Oceanic explicit time stepping
    ########################################################################

    if cpl==1

        TOst[1]=TO[i-1,1]+(SW[i-1,1]-LH[i-1]-SH[i-1]-LW[i-1])/cp/(dz*rhow[1])*dt-
            (1-Oimp)*(kot[1]*(TO[i-1,1]-TO[i-1,2])/dz^2*dt);
        TOst[end]=TO[i-1,end]+SW[i-1,end]/cp/(dz*rhow[end])*dt+
            (1-Oimp)*(kot[end-1].*(TO[i-1,end-1]-TO[i-1,end])/dz^2*dt);
        TOst[2:end-1]=TO[i-1,2:end-1]+(1-Oimp)*(
            kot[1:end-2].*(TO[i-1,1:end-2]-TO[i-1,2:end-1])/dz^2*dt-
            kot[2:end-1].*(TO[i-1,2:end-1]-TO[i-1,3:end  ])/dz^2*dt)+
            SW[i-1,2:end-1]./cp./(dz*rhow[2:end-1])*dt;

        UOst[1]=UO[i-1,1]+
             (TAU[i-1])*TAUS/(rhow[1])/dz*dt-
             (1-Oimp)*(kom[1]*(UO[i-1,1]-UO[i-1,2])/dz^2*dt);
        UOst[end]=UO[i-1,end]+
             (1-Oimp)*(kom[end-1].*(UO[i-1,end-1]-UO[i-1,end])/dz^2*dt);
        UOst[2:end-1]=UO[i-1,2:end-1]+(1-Oimp)*(
             kom[1:end-2].*(UO[i-1,1:end-2]-UO[i-1,2:end-1])/dz^2*dt-
             kom[2:end-1].*(UO[i-1,2:end-1]-UO[i-1,3:end  ])/dz^2*dt);

    end

    if Aimp==0
        UA[i,:]=UAst; qA[i,:]=qAst; THETA[i,:]=THETAst;
    else
        du=-Aimp*kam[1:end-1].*dt./H^2; dl=-Aimp*kam[2:end].*dt./H^2;
        d=(1 .+Aimp.*dt./H^2 .*([kam[1:end-1]; 0]+[0; kam[2:end]]));
        A=inv(Tridiagonal(dl,d,du));
        UA[i,:]=A*UAst;
        du=-Aimp*kat[1:end-1].*dt./H^2; dl=-Aimp*kat[2:end].*dt./H^2;
        d=(1 .+Aimp.*dt./H^2 .*([kat[1:end-1]; 0]+[0; kat[2:end]]));
        A=inv(Tridiagonal(dl,d,du));
        qA[i,:]=A*qAst; THETA[i,:]=A*THETAst;
    end

    if Oimp==0
        UO[i,:]=UOst; TO[i,:]=TOst;
    else
        du=-Oimp*kom[1:end-1].*dt./dz^2; dl=-Oimp*kom[2:end].*dt./dz^2;
        d=(1 .+Oimp.*dt./dz^2 .*([kom[1:end-1]; 0]+[0; kom[2:end]]));
        A=inv(Tridiagonal(dl,d,du));
        UO[i,:]=A*UOst;
        du=-Oimp*kot[1:end-1].*dt./dz^2; dl=-Oimp*kot[2:end].*dt./dz^2;
        d=(1 .+Oimp.*dt./dz^2 .*([kot[1:end-1]; 0]+[0; kot[2:end]]));
        A=inv(Tridiagonal(dl,d,du));
        TO[i,:]=A*TOst;
    end
    SST[i]=TO[i,1]

    TA[i,:]=(THETA[i,:].+d2k).*(P./P[1]).^(2/7).-d2k;
    rhoa[i,:]=P./(Rgas.*(TA[i,:].+d2k));

    if moist==1
        qsat=saltsat*cvapor_fac*exp.(-cvapor_exp./(TA[i,:].+d2k))./rhoa[i-1,:];
        dq=(qA[i,:].-qsat).*((qA[i,:].-qsat).>0);
        dq=min.(dq,qA[i,:]);
        dTHETA=Av.*dq./cp; #[J/kg][j-1 kg K]
        #println([i maximum(dTHETA) maximum(dq)])
        ind=qA[i,:].>qsat;
        qA[i,ind]=qA[i,ind]-dq[ind];
        #THETA[i,:]=THETA[i,: ]+dTHETA;
    end

    if rem(i,3600*6/dt)==0;#3600*6/dt
      println(i)
      f1=plot(-ZO[1:end-1],TO[i,:],ylabel=("TO"));#,ylims=(19,21)
      f2=plot(-ZO[1:end-1],UO[i,:],ylabel=("UO"));
      f3=plot(-ZO[1:end-1],kom,ylabel=("KOm"));
      f4=plot(ZA,UA[i,:],ylabel=("UA"));
      f5=plot(ZA,Tv.-273.15,ylabel=("Tv")); #ylim([19 21]);...
      f6=plot(ZA,qA[i,:],ylabel=("q"));
      f7=plot(ZA,kam,ylabel=("KAm"));
      plt=plot(f1,f2,f3,f4,f5,f6,f7,layout=(7,1),legend=false,titlefontsize=6)
      display(plt)
    end

end

k1=plot(rollmean(U0[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="U0");
k2=plot(rollmean(UA[:,end],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="Utop");
k3=plot(rollmean(UA[:,end]-U0[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="Utop-U0");
k4=plot(rollmean(UA[:,1],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="U");
k5=plot(rollmean(THETA[:,2],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="THETA");
k6=plot(rollmean(SST,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="SST");
k7=plot(rollmean(TO[:,10],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              ylabel="T");
#tt=string("DZ=",dz)
plt=plot(k1,k2,k3,k4,k5,k6,k7,layout=(7,1),legend=false)
display(plt)

#PyPlot.suptitle(tt)

k1=plot(rollmean(xi[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="xi");
k2=plot(rollmean(ustar[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="ustar");
k3=plot(rollmean(psimh[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="psimh");
k4=plot(rollmean(psixh[1,:],Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="psixh");
k5=plot(rollmean((rd[1,:].^2),Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="CD");
k6=plot(rollmean((rh[1,:].*rd[1,:]),Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="CH");
k7=plot(rollmean((re[1,:].*rd[1,:]),Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="CE");
#tt=string("DZ=",dz)
plt=plot(k1,k2,k3,k4,k5,k6,k7,layout=(7,1),legend=false)
display(plt)

#PyPlot.suptitle(tt)

k1=plot(rollmean(SW[:,1]/cp/(dz*rhow0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="SW");
k2=plot(rollmean(LW[1,:]/cp/(dz*rhow0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="LW");
k3=plot(rollmean(SH[1,:]/cp/(dz*rhow0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="SH");
k4=plot(rollmean(LH[1,:]/cp/(dz*rhow0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="LH");
k5=plot(rollmean((SW[:,1]-LW[1,:]-SH[1,:]-LH[1,:])/cp/(dz*rhow0)*dt,Int(24*60*60/dt)),
              xticks=(0:24*60*60/dt:n,string.(0:ndays)),
              title="QNET");
plt=plot(k1,k2,k3,k4,k5,layout=(5,1),legend=false);
display(plt)

#PyPlot.suptitle(tt)
