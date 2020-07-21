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

host=NFW('mv',mv0,'cc',cv0,'fg',fg0);

R0=1.5.*host.Rvir;
rp=logspace(-2,0,1000).*R0;
rhoICM=gasDensity(host,rp,'kpc');
rp0=rp;        


%% prepare changeing parameters

len=30;
cmap=brewermap(len,'Spectral');
ticInd=round(linspace(1,len,10));


vc=vcirc(host,R0,'kpc');

vxP=host.Vvir.*linspace(0,1,len);

vyP=vc.*linspace(0,1,len);


str={'vx','vy'}; %,'cv','fg'};


%% do a fiducal run


% do orbit

radIC.x=1.5.*host.Rvir;
radIC.y=0;
radIC.vx=0;%-host.Vvir;%./sqrt(2);
radIC.vy=0;%host.Vvir./sqrt(2);

tOrbit=2*pi*host.Rvir./host.Vvir;
tmax=4*tOrbit;

dt=tOrbit/1e4;
dtmin=dt/1e3;

radOrb=orbits.rk4_orbitIntegration(radIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
timeunit=(Units.kpc/Units.km)/(Units.Gyr);


% set pram

%v2=radOrb.vx.^2+radOrb.vy.^2;

%dr=diff(rp);
ind=find(diff(radOrb.rad)>0,1,'first')-1;

vsat=interp1(radOrb.rad(1:ind),radOrb.vel(1:ind),rp,'pchip');

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
for j=2:length(str)
    
    hf=figure;
    
    
    mstrip=zeros(length(rp));
    
    
    radIC.vx=0;
    radIC.vy=0;
    
    for i=1:len
        
        
        switch str{j}
            case 'vx'
                radIC.vx=-vxP(i);
            case 'vy'
                radIC.vy=vyP(i);
        end
        
        
        %% define host
        
        tOrbit=2*pi*(host.Rvir)./host.Vvir;
        tmax=4*tOrbit;
        
        dt=tOrbit/1e4;
        dtmin=dt/1e3;
        
        radOrb=orbits.rk4_orbitIntegration(radIC,dt,tmax,dtmin,@orbits.rhs_nfw,host);
                
        
        % set pram
        
        %v2=radOrb.vx.^2+radOrb.vy.^2;
        
        %dr=diff(rp);
        ind=find(diff(radOrb.rad)>0,1,'first');
        if isempty(ind)
            ind=length(radOrb.rad);
        elseif ind==1
                    ind=find(diff(radOrb.rad)>0,2,'first');
        ind=ind(2);
        end
        
        
        rp=logspace(log10(radOrb.rad(ind)),log10(R0),1000);
        rhoICM=gasDensity(host,rp,'kpc');

        vsat=interp1(radOrb.rad(1:ind),radOrb.vel(1:ind),rp,'pchip');
        
        pram=alf0.*rhoICM.*vsat.^2;
        
        %% find stripping radius and stripped mass profile
        
        rstrip=zeros(size(rp));
        
        for k=length(rp):-1:1
            
            
            
            
            pr=pram(k)/fmax;
            
            if pr<1
                rstrip(k)=interp1(ff(fmaxInd:end)./fmax,rg(fmaxInd:end),pr,'pchip');
                
            end
            
        end
        
                
        %mstrip(i,:)=mass(gal.GasDisk,rstrip,'kpc');
        mstrip=mass(gal.GasDisk,rstrip,'kpc');
        
        fprintf('%s %% finished \n',num2str(i/len*100));
        
        plot(rp./host.Rvir,mstrip./gal.GasDisk.Md,'color',cmap(i,:))
        scatter(rp(1)./host.Rvir,mstrip(1)./gal.GasDisk.Md,40,cmap(i,:),'filled')
        
        
        if i==1
            hold on
        end
        
        
        
    end
    
    
    % plot fiducial case
    plot(rp0./host.Rvir,mstrip1(:)./gal.GasDisk.Md,'color','k','linewidth',1.5)
    
    %% plot figure
    switch str{j}
        case 'vx'
            tic=vxP(ticInd)./host.Vvir;
            barLab='$V_x/V_\mathrm{vir}$';
        case 'vy'
            tic=vyP(ticInd)./vc;
            barLab='$V_y/V_\mathrm{circ}$';
        
    end
    
    for k=1:length(tic)
        lab{k}=num2str(tic(k));
    end
    
    grid
    xlim([0.01 1.52])
    ylim([0 1.02])
    
    set(gca,'color',[0.9 0.9 0.9])
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
