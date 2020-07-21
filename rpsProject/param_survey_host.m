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

alf0=0.5;

rg=(0.01:0.01:10).*gal.GasDisk.Rd;
gasForce=gasBindingForce(gal,rg,'kpc');

ff=abs(gasForce);
[fmax,fmaxInd]=max(ff);


%% prepare host set parameters

mv0=1e14;
cv0=5;
fg0=0.15;


%% prepare changeing parameters

len=30;
cmap=brewermap(len,'Spectral');
ticInd=round(linspace(1,len,10));

mvRange=[13.5 15.5];

mvP=10.^linspace(mvRange(1),mvRange(2),len);

cvRange=[2 15];
cvP=logspace(log10(cvRange(1)),log10(cvRange(2)),len);

fgRange=[0.05 0.25];
fgP=logspace(log10(fgRange(1)),log10(fgRange(2)),len);


str={'mv','cv','fg'};



%% do a fiducal run

% set host

host=NFW('mv',mv0,'cc',cv0,'fg',fg0);

% do orbit

radIC.x=1.5.*host.Rvir;
radIC.y=0;
radIC.vx=0; %-host.Vvir;%./sqrt(2);
radIC.vy=0;%host.Vvir./sqrt(2);

tOrbit=2*pi*host.Rvir./host.Vvir;
tmax=3*tOrbit;

dt=tOrbit/1e4;
dtmin=dt/1e3;

radOrb=orbits.rk4_orbitIntegration(radIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
timeunit=(Units.kpc/Units.km)/(Units.Gyr);


% set pram
rp=logspace(-2,0,1000).*radOrb.rad(1);

v2=radOrb.vx.^2+radOrb.vy.^2;

%dr=diff(rp);
ind=find(diff(radOrb.rad)>0,1,'first')-1;

vsat=interp1(radOrb.rad(1:ind),radOrb.vel(1:ind),rp,'pchip');
rhoICM=gasDensity(host,rp,'kpc');
pram=alf0.*rhoICM.*vsat.^2;

%% find stripping radius and stripped mass profile

rstrip=zeros(size(rp));

for k=length(rp):-1:1
    
    pr=pram(k)/fmax;
    
    if pr<1
        rstrip(k)=interp1(ff(fmaxInd:end)./fmax,rg(fmaxInd:end),pr,'pchip');
        
    end
    
end

mstrip1=mass(gal.GasDisk,rstrip,'kpc');

%% do the parameter run
for j=1:length(str)
    
    hf=figure;
    
    
    mstrip=zeros(len,length(rp));
    
    mv=mv0;cv=cv0;fg=fg0;
    
    
    for i=1:len
        
        
        switch str{j}
            case 'mv'
                mv=mvP(i);
            case 'fg'
                fg=fgP(i);
            case 'cv'
                cv=cvP(i);
        end
        
        
        %% define host
        
        % set host
        
        host=NFW('mv',mv,'cc',cv,'fg',fg);
        
        % do orbit
        
        radIC.x=1.5.*host.Rvir;
        
        tOrbit=2*pi*(host.Rvir)./host.Vvir;
        tmax=3*tOrbit;
        
        dt=tOrbit/1e4;
        dtmin=dt/1e3;
        
        radOrb=orbits.rk4_orbitIntegration(radIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
                
        
        % set pram
        rp=logspace(-2,0,1000).*radOrb.rad(1);
        
        v2=radOrb.vx.^2+radOrb.vy.^2;
        
        %dr=diff(rp);
        ind=find(diff(radOrb.rad)>0,1,'first')-1;
        
        vsat=interp1(radOrb.rad(1:ind),radOrb.vel(1:ind),rp,'pchip');
        rhoICM=gasDensity(host,rp,'kpc');
        pram=alf0.*rhoICM.*vsat.^2;
        
        %% find stripping radius and stripped mass profile
        
        rstrip=zeros(size(rp));
        
        for k=length(rp):-1:1
            
            pr=pram(k)/fmax;
            
            if pr<1
                rstrip(k)=interp1(ff(fmaxInd:end)./fmax,rg(fmaxInd:end),pr,'pchip');
                
            end
            
        end
        
                
        mstrip(i,:)=mass(gal.GasDisk,rstrip,'kpc');
        
        fprintf('%s %% finished \n',num2str(i/len*100));
        
        plot(rp./host.Rvir,mstrip(i,:)./gal.GasDisk.Md,'color',cmap(i,:))
        
        
        if i==1
            hold on
        end
        
        
        
    end
    
    
    % plot fiducial case
    plot(rp./host.Rvir,mstrip1(:)./gal.GasDisk.Md,'color','k','linewidth',1.5)
    
    %% plot figure
    switch str{j}
        case 'mv'
            tic=mvP(ticInd);
            barLab='$M_\mathrm{host}$';
        case 'cv'
            tic=cvP(ticInd);
            barLab='$c_\mathrm{vir,H}$';
        case 'fg'
            tic=fgP(ticInd);
            barLab='$f_\mathrm{gas}$';
    end
    
    for k=1:length(tic)
        lab{k}=num2str(tic(k));
    end
    
    grid
    xlim([0.01 1.52])
    ylim([0 1.02])
    
    %set(gca,'color','k')
    xlabelmine('$r_\mathrm{pos}/R_\mathrm{vir}$')
    ylabelmine('$M(<r_\mathrm{rstrip})/M_{\mathrm{gas}}$')
    
    hb=colorbar('colormap',cmap,'ticks',linspace(0,1,10),'ticklabels',lab);
    barTitle(hb,barLab);
    
    
    %
    %
    %     %% plotting
    %     figure
    %
    %     h(1)=plot(rppR./host.Rvir,mstripRad./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Radial');
    %     hold on
    %     h(2)=plot(rppEcc./host.Rvir,mstripEcc./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Eccentric');
    %
    %
    %
    %     hl=legend(h);
    %     set(hl,'interpreter','latex','location','SouthEast','fontsize',14)
    %
    %     xlim([0 1.1])
    %     ylim([0 1])
    %
    %
    %
    %     % vs time
    %     figure
    %
    %
    %
    %     h(1)=plot(orbRad.t(1:indR).*timeunit,mstripRad./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Radial');
    %     hold on
    %     h(2)=plot(orbEcc.t(1:indE).*timeunit,mstripEcc./gal.GasDisk.Md,'linewidth',1.5,'DisplayName','Eccentric');
    %
    %     grid
    %
    %     hl=legend(h);
    %     set(hl,'interpreter','latex','location','SouthEast','fontsize',14)
    %
    %     xlim([0 1.21])
    %     ylim([0 1])
    %
    %
    %     xlabelmine('$t\,[\mathrm{Gyr}]$')
    %     ylabelmine('$M(<r_\mathrm{rstrip})/M_{\mathrm{gas}}$')
    %
    %
    %
end
