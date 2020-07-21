%% We carry out a parameter survey of galaxy parameters without changeing the total mass: beta, xi and


    % perliminaries
    units;
    
    %% prepare sat aarametes
    fgs=0.25;
    fbs=0.15;
    
    Ms=1e10;
    rd=5;
    %Mg=fgs.*Ms;
    beta=1;
    %Mb=fbs.*Ms;
    xi=3;
    Mh=1e11;
    cs=5;
    
    
    gal=GALAXY('ms',Ms,'rd',rd,'fgs',fgs,'beta',beta,'fbs',fbs,'xi',xi,'Mh',Mh,'cv',cs);
    gal2=GALAXY('ms',Ms,'rd',rd,'fgs',fgs,'beta',beta,'fbs',fbs,'xi',xi,'Mh',Mh,'cv',cs);
    
    % alf0=0.5;
    %
    % rg=(0.01:0.01:10).*gal.GasDisk.Rd;
    % gasForce=gasBindingForce(gal,rg,'kpc');
    %
    % ff=abs(gasForce);
    % [fmax,fmaxInd]=max(ff);
    
    
    %% prepare host set parameters
    
    mv0=1e14;
    cv0=5;
    fg0=0.15;
    
     % set host
    
    host=NFW('mv',mv0,'cc',cv0,'fg',fg0);
    


%if orbitFlag
    %% do a fiducal run
    
   
    
    % do orbit
    
    radIC.x=1.5.*host.Rvir;
    radIC.y=0;
    radIC.vx=0; %-host.Vvir;%./sqrt(2);
    radIC.vy=0; %host.Vvir./sqrt(2);
    
    orbTag='rad';
    
    
    tOrbit=2*pi*host.Rvir./host.Vvir;
    tmax=3*tOrbit;
    
    dt=tOrbit/1e4;
    dtmin=dt/1e3;
    
    radOrb=orbits.rk4_orbitIntegration(radIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
    timeunit=(Units.kpc/Units.km)/(Units.Gyr);
%end

%% run over orbit without rps
%if sfFlag
    len=length(radOrb.t);
    gmass=zeros(length(gal.gasMass.mass),len);
    mgas=zeros(1,len);
    mstar=zeros(1,len);
    sfr1=zeros(length(gal.gasMass.mass),len);
    
    
    gmass(:,1)=gal.gasMass.mass;
    mstar(1)=gal.StellarDisk.Md;
    mgas(1)=gal.GasDisk.Md;
    sfr1(:,1)=GALAXY.starFormationRate(gal.gasMass.rr,gal.gasMass.mass);
    
    for i=2:length(radOrb.t)
        
        
        dt=diff(radOrb.t(i-1:i)).*timeunit; % in Gyr
        
        % find out how much gas is converted to stars
        sfr=GALAXY.starFormationRate(gal.gasMass.rr,gal.gasMass.mass);
        sfr1(:,i)=sfr;
        dmass=sfr.*(dt.*1e9);
        
        mstar(i)=mstar(i-1)+sum(dmass);
        
        % update the gas mass.
        gal.gasMass.mass=max(gal.gasMass.mass-dmass,zeros(size(dmass)));
        gmass(:,i)=gal.gasMass.mass;
        
        mgas(i)=sum(gmass(:,i));
    end
    
%end
%% run again with RPS

%if rpFlag
    alf=1;
    ind=find(diff(radOrb.rad)>0,1,'first')-1;
    rp=radOrb.rad(1:ind);
    
    
    rhoICM=gasDensity(host,rp,'kpc');
    vsat2=radOrb.vx(1:ind).^2+radOrb.vy(1:ind).^2;
    pram0=alf.*rhoICM.*vsat2;
    
    
    rg=(0.01:0.01:10).*gal2.GasDisk.Rd;
    gasForce=gasBindingForce(gal2,rg,'kpc');
    
    ff=abs(gasForce);
    [fmax,fmaxInd]=max(ff);
    
    
    
    len=ind;
    gmass2=zeros(length(gal2.gasMass.mass),len);
    mgas2=zeros(1,len);
    mstar2=zeros(1,len);
    sfr2=zeros(length(gal2.gasMass.mass),len);
    
     pr=pram0(1)/fmax;
     gmass2(:,1)=gal2.gasMass.mass;
     mstar2(1)=gal2.StellarDisk.Md;
     mgas2(1)=gal2.GasDisk.Md;
     sfr2(:,1)=GALAXY.starFormationRate(gal2.gasMass.rr,gal2.gasMass.mass);
     rstrip=zeros(size(rp));
     rstrip(1)=interp1(ff(fmaxInd:end)./fmax,rg(fmaxInd:end),pr,'pchip');
         
    for i=2:length(rp)
        
        %     ff=abs(gasForce);
        %     [fmax,fmaxInd]=max(ff);
        %
        
        pr=pram0(i)/fmax;
        
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
        
        stripInd=find(gal2.gasMass.rr(2:end)>rstrip(i));
        gal2.gasMass.mass(stripInd)=0;
        
        
        dt=diff(radOrb.t(i-1:i)).*timeunit; % in Gyr
        
        % find out how much gas is converted to stars
        sfr=GALAXY.starFormationRate(gal2.gasMass.rr,gal2.gasMass.mass);
        sfr2(:,i)=sfr;
        dmass=sfr.*(dt.*1e9);
        
        mstar2(i)=mstar2(i-1)+sum(dmass);
        
        % update the gas mass.
        gal2.gasMass.mass=max(gal2.gasMass.mass-dmass,zeros(size(dmass)));
        gmass2(:,i)=gal2.gasMass.mass;
        
        mgas2(i)=sum(gmass2(:,i));
        
    end
    
    
    
%end

%% figures 
tim=radOrb.t(1:ind).*timeunit;
timLab='Time [Gyr]'; 
% plot mstar & mgas 

cmap=brewermap(8,'Set1');

h=[];
figure 
h(1)=plot(tim,mgas(1:ind)./mgas(1),'color',cmap(1,:),'DisplayName','SF Only');
hold on 
h(2)=plot(tim,mgas2./mgas2(1),'color',cmap(2,:),'DisplayName','SF \& RPS');

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)


grid
xlabelmine(timLab);
ylabelmine('Gas Mass Loss');%$[\mathrm{M_\odot}]$')
set(gca,'Fontsize',14)

printout_fig(gcf,['sfRps_' orbTag '_compare_mgas'],'v')

% ------------------------------

h=[];
figure 
h(1)=plot(tim,mstar(1:ind)./mstar(1),'color',cmap(1,:),'DisplayName','SF Only');
hold on 
h(2)=plot(tim,mstar2./mstar2(1),'color',cmap(2,:),'DisplayName','SF \& RPS');

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14,'location','SouthEast')


grid
xlabelmine(timLab);
ylabelmine('Stellar Mass Growth'); %$[\mathrm{M_\odot}]$')
set(gca,'Fontsize',14)

printout_fig(gcf,['sfRps_' orbTag '_compare_mstar'],'v')

% ------------------------------

h=[];
figure 
h(1)=plot(tim,sum(sfr1(:,1:ind),1)./sum(sfr1(:,1),1),'color',cmap(1,:),'DisplayName','SF Only');
hold on 
h(2)=plot(tim,sum(sfr2,1)./sum(sfr2(:,1),1),'color',cmap(2,:),'DisplayName','SF \& RPS');

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)


grid
xlabelmine(timLab);
ylabelmine('SFR Decay'); %$[\mathrm{M_\odot\,yr^{-1}}]$')
set(gca,'Fontsize',14)

printout_fig(gcf,['sfRps_' orbTag '_compare_sfr'],'v')

% ------------------------------

h=[];
figure 
h(1)=plot(tim,sum(sfr1(:,1:ind),1)./mstar(1:ind),'color',cmap(1,:),'DisplayName','SF Only');
hold on 
h(2)=plot(tim,sum(sfr2,1)./mstar2,'color',cmap(2,:),'DisplayName','SF \& RPS');
plot(tim,ones(size(tim)).*1e-11,'--k')

hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)


grid
xlabelmine(timLab);
ylabelmine('sSFR $[\mathrm{yr^{-1}}]$')
set(gca,'Fontsize',14)

printout_fig(gcf,['sfRps_' orbTag '_compare_ssfr'],'v')

% ------------------------------

db=floor(ind/12);
ii=1:db:ind;
ii(end)=ind;

figure
h=[];
cnt=0;
for i=ii
    cnt=cnt+1;
    tag=sprintf('%3.2f Gyr',tim(i));
    h(cnt)=plot(gal.gasMass.rr(2:end),gmass(:,i),'Displayname',tag);
    hold on
end
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)

grid
xlabelmine('galactic radius [kpc]');
ylabelmine('Gas Mass $[\mathrm{M_\odot}]$')
set(gca,'Fontsize',14)
titlemine('SF Only')

printout_fig(gcf,['sf_' orbTag '_gasProf'],'v')

figure
h=[];
cnt=0;
for i=ii
        cnt=cnt+1;
    tag=sprintf('%3.2f Gyr',tim(i));
    h(cnt)=plot(gal2.gasMass.rr(2:end),gmass2(:,i),'Displayname',tag);
    hold on
end
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)

grid
xlabelmine('galactic radius [kpc]');
ylabelmine('Gas Mass $[\mathrm{M_\odot}]$')
set(gca,'Fontsize',14)
titlemine('SF \& RPS')

printout_fig(gcf,['sfRPS_' orbTag '_gasProf'],'v')


% -----------------------


    figure
h=[];
cnt=0;
for i=ii
    cnt=cnt+1;
    tag=sprintf('%3.2f Gyr',tim(i));
    h(cnt)=plot(gal.gasMass.rr(2:end),sfr1(:,i),'Displayname',tag);
    hold on
end
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)

grid
xlabelmine('galactic radius [kpc]');
ylabelmine('SFR $[\mathrm{M_\odot\,yr^{-1}}]$')
set(gca,'Fontsize',14)
titlemine('SF Only')

printout_fig(gcf,['sf_' orbTag '_sfrProf'],'v')

figure
h=[];
cnt=0;
for i=ii
        cnt=cnt+1;
    tag=sprintf('%3.2f Gyr',tim(i));
    h(cnt)=plot(gal2.gasMass.rr(2:end),sfr2(:,i),'Displayname',tag);
    hold on
end
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14)

grid
xlabelmine('galactic radius [kpc]');
ylabelmine('SFR $[\mathrm{M_\odot \,yr^{-1}}]$')
set(gca,'Fontsize',14)
titlemine('SF \& RPS')

printout_fig(gcf,['sfRPS_' orbTag '_sfrProf'],'v')


% -----------------------
    




% plot gas mass profile 











