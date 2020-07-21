%% We carry out a parameter survey of galaxy parameters without changeing the total mass: beta, xi and

% perliminaries
units;

%% set host
host=NFW('mv',1e14,'cc',5,'fg',0.15);

%% prepare orbits
if doOrbit
    radIC.x=1.5.*host.Rvir;
    radIC.y=0;
    radIC.vx=0; %-host.Vvir;%./sqrt(2);
    radIC.vy=0;%host.Vvir./sqrt(2);
    
    % eccIC.x=1.*host.Rvir;
    % eccIC.y=0;
    % eccIC.vx=-host.Vvir./sqrt(2);
    % eccIC.vy=host.Vvir./sqrt(2);
    
    
    tOrbit=2*pi*sqrt(host.Rvir)./host.Vvir;
    tmax=10*tOrbit;
    
    dt=tOrbit/1e5;
    dtmin=dt/1e3;
    
    radOrb=orbits.rk4_orbitIntegration(radIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
    % eccOrb=orbits.rk4_orbitIntegration(eccIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
    
    timeunit=(Units.kpc/Units.km)/(Units.Gyr);
    %time=orb.t.*timeunit;
    
end
%% prepare set parametes
fgs=0.25;
fbs=0.15;

Ms=1e10;
rd=5;
Mg=fgs.*Ms;
beta0=1;
Mb=fbs.*Ms;
xi0=3;
Mh=1e11;
cv0=5;

alf0=0.5;

%% prepare changeing parameters

len=30;
cmap=brewermap(len,'Spectral');
ticInd=round(linspace(1,len,10));

betaRange=[0.2 5];

betaP=logspace(log10(betaRange(1)),log10(betaRange(2)),len);

xiRange=[0.667 5];
xiP=logspace(log10(xiRange(1)),log10(xiRange(2)),len);

cvRange=[2 15];
cvP=logspace(log10(cvRange(1)),log10(cvRange(2)),len);

alfRange=[0.05 2];
alfP=logspace(log10(alfRange(1)),log10(alfRange(2)),len);


%rp=(0.01:0.01:1).*max(radOrb.rad);
rp=logspace(-2,0,1000).*max(radOrb.rad);

v2=radOrb.vx.^2+radOrb.vy.^2;

%dr=diff(rp);
ind=find(diff(radOrb.rad)>0,1,'first')-1;

vsat=interp1(radOrb.rad(1:ind),radOrb.vel(1:ind),rp,'pchip');
rhoICM=gasDensity(host,rp,'kpc');
pram0=rhoICM.*vsat.^2;

str={'beta','xi','cv','alf'};


%% do one run with fiducial values
%% define galaxy
gal=GALAXY('ms',Ms,'rd',rd,'fgs',fgs,'beta',beta0,'fbs',fbs,'xi',xi0,'Mh',Mh,'cv',cv0);

%% calculate the binding force
rg=(0.01:0.01:10).*gal.GasDisk.Rd;
gasForce=gasBindingForce(gal,rg,'kpc');

ff=abs(gasForce);
[fmax,fmaxInd]=max(ff);
%% calculate RPS

pram=alf0.*pram0;

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
    for i=1:len
        
        beta=beta0;xi=xi0;cv=cv0;alf=alf0;
        
        switch str{j}
            case 'beta'
                beta=betaP(i);
            case 'xi'
                xi=xiP(i);
            case 'cv'
                cv=cvP(i);
            case 'alf'
                alf=alfP(i);
        end
        
        
        %% define galaxy
        gal=GALAXY('ms',Ms,'rd',rd,'fgs',fgs,'beta',beta,'fbs',fbs,'xi',xi,'Mh',Mh,'cv',cv);
        
        %% calculate the binding force
        rg=(0.01:0.01:10).*gal.GasDisk.Rd;
        gasForce=gasBindingForce(gal,rg,'kpc');
        
        ff=abs(gasForce);
        [fmax,fmaxInd]=max(ff);
        %% calculate RPS
        
        
        
        
        pram=alf.*pram0;
        
        %% find stripping radius and stripped mass profile
        
        rstrip=zeros(size(rp));
        
        
        for k=length(rp):-1:1
            
            
            
            
            pr=pram(k)/fmax;
            
            if pr<1
                rstrip(k)=interp1(ff(fmaxInd:end)./fmax,rg(fmaxInd:end),pr,'pchip');
%                 %                 if i>1
%                 %                     rstrip(k)=min(rs,rstrip(k-1));
%                 %                 else
%                 rstrip(k)=rs;
%                 %                 end
%             else
%                 rstrip(k)=0;
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
        case 'beta'
            tic=betaP(ticInd);
            barLab='$R_\mathrm{s}/R_\mathrm{g}$';
        case 'xi'
            xi=xiP(i);
             barLab='$R_\mathrm{s}/R_\mathrm{b}$';
        case 'cv'
            tic=cvP(ticInd);
             barLab='$c_\mathrm{vir}$';
        case 'alf'
            tic=alfP(ticInd);
             barLab='$\epsilon$';
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
