%% define stuff

units;

host=NFW('mv',1e14,'cc',5,'fg',0.15);

gal=GALAXY('ms',1e10,'rd',5,'fgs',0.25,'beta',1,'fbs',0.05,'xi',1,'Mh',1e11,'cv',5);


alf=0.5;

%% prepare orbit 
IC.x=1.*host.Rvir;
IC.y=0;
IC.vx=-host.Vvir;%./sqrt(2);
IC.vy=0;%host.Vvir./sqrt(2);

tOrbit=2*pi*sqrt(IC.x.^2+IC.y^2)./sqrt(IC.vx.^2+IC.vy.^2);
tmax=3*tOrbit;

dt=tOrbit/1e4;
dtmin=dt/1e3;

orb=orbits.rk4_orbitIntegration(IC,dt,tmax,dtmin,@orbits.rhs_nfw,host);

timeunit=(Units.kpc/Units.km)/(Units.Gyr);
time=orb.t.*timeunit;


%% calculate RPS 

rp=orb.rad;
dr=diff(rp);
ind=find(dr>0,1,'first')-1;


rhoICM=gasDensity(host,rp,'kpc');

vsat2=orb.vx.^2+orb.vy.^2;

pram0=alf.*rhoICM.*vsat2;

%% calculate the binding force 
rg=(0.01:0.01:10).*gal.GasDisk.Rd;

gasForce=gasBindingForce(gal,rg,'kpc');

%% find stripping radius and stripped mass profile 

rpp=rp(1:ind);
pram=pram0(1:ind);
rstrip=zeros(size(rpp)

for i=1:length(rpp)

    ff=abs(gasForce);
    [fmax,fmaxInd]=max(ff);
    
        
    pr=pram(i)./fmax;
    
    if pr<1
        
        rs=interp1(ff(fmaxInd:end)./fmax,rg(fmaxInd:end),pr,'pchip');
        if i>1
            rstrip(i)=min(rs,rstrip(i-1));
        else
            rstrip(i)=rs;
        end
    else
        rstrip(i)=0;
    end
    
    
    %    
%     ff=abs(gasForce)-pram(i);
%     [fmax,fmaxInd]=max(ff);
%     
%     ff(1:fmaxInd)=fmax;
%     
%     rs=interp1(ff,rg,0,'pchip');
%    if i>1
%     rstrip(i)=min(rs,rstrip(i-1));
%    else
%        rstrip(i)=rs;
%    end
%    
%    
   
end
    
mstrip=mass(gal.GasDisk,rstrip,'kpc');


%% plotting 
figure

h(1)=plot(rppR./host.Rvir,mstripRad./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Radial');
hold on
h(2)=plot(rppEcc./host.Rvir,mstripEcc./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Eccentric');

grid

hl=legend(h);
set(hl,'interpreter','latex','location','SouthEast','fontsize',14)

xlim([0 1.1])
ylim([0 1])

xlabelmine('$r_\mathrm{pos}/R_\mathrm{vir}$')
ylabelmine('$M(<r_\mathrm{rstrip})/M_{\mathrm{gas}}$')


% vs time 
figure



h(1)=plot(orbRad.t(1:indR).*timeunit,mstripRad./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Radial');
hold on
h(2)=plot(orbEcc.t(1:indE).*timeunit,mstripEcc./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Eccentric');

grid

hl=legend(h);
set(hl,'interpreter','latex','location','SouthEast','fontsize',14)

xlim([0 1.21])
ylim([0 1])


xlabelmine('$t\,[\mathrm{Gyr}]$')
ylabelmine('$M(<r_\mathrm{rstrip})/M_{\mathrm{gas}}$')

