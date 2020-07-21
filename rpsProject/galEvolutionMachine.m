function res=galEvolutionMachine(gal,host,orb,varargin)

units;
timeunit=(Units.kpc/Units.km)/(Units.Gyr);
time=orb.t.*timeunit;
% defaults
alf=1;

wbarFlag=true;

sfrAlfa=1.4;
sfrAfac=2.5;
plotFlag=false;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'alf','alfa','rpsalpha'}
            i=i+1;
            alf=varargin{i};
        case{'sfra','sfrafac'}
            i=i+1;
            sfrAfac=varargin{i};
        case{'sfralfa','sfralf'}
            i=i+1;
            sfrAlfa=varargin{i};
        case 'plot'
            plotFlag=true;
        case {'nobar'}
            wbarFlag=false;
        otherwise
            error('GALEVOLUTIONMACHINE - Illegal argument: %s',varargin{i});
    end
    i=i+1;
    
end



% prepare host side of RPS
rp=orb.rad;
rhoICM=gasDensity(host,rp,'kpc');
vsat2=orb.vx.^2+orb.vy.^2;
pram0=alf.*rhoICM.*vsat2;

% prepare galaxy side of RPS
rgal=(0:0.01:50).*gal.GasDisk.Rd; % in kpc
gasForce=gasBindingForce(gal,rgal,'kpc');
ff=abs(gasForce);
[fmax,fmaxInd]=max(ff);


%initialize
len=length(rp);
gasProf=zeros(length(gal.gasMass.mass),len);
gasMass=zeros(1,len);
stellarMass=zeros(1,len);
sfrProf=zeros(length(gal.gasMass.mass),len);

% do first index
pr=pram0(1)/fmax;
gasProf(:,1)=gal.gasMass.mass;
stellarMass(1)=gal.StellarDisk.Md;
gasMass(1)=gal.GasDisk.Md;
sfrProf(:,1)=GALAXY.starFormationRate(gal.gasMass.rr,gal.gasMass.mass,...
    'alf',sfrAlfa,'afac',sfrAfac);

rstrip=zeros(size(rp));
rstrip(1)=interp1(ff(fmaxInd:end)./fmax,rgal(fmaxInd:end),pr,'pchip');


if wbarFlag; hb = waitbar(0,'Galaxy is evolving...');end

for i=2:len
    
    if wbarFlag
        if mod(i,100)==0
            waitbar(i./len,hb)
        end
    end
    
    if alf>0
    pr=pram0(i)/fmax;
    
    if pr<1
        
        rs=interp1(ff(fmaxInd:end)./fmax,rgal(fmaxInd:end),pr,'pchip');
        rstrip(i)=min(rs,rstrip(i-1));
        
    else
        rstrip(i)=0;
    end
    
    stripMask=gal.gasMass.rr(2:end)>rstrip(i);
    gal.gasMass.mass(stripMask)=0;
    end
    dt=diff(time(i-1:i)); % in Gyr
    
    % find out how much gas is converted to stars
    sfr=GALAXY.starFormationRate(gal.gasMass.rr,gal.gasMass.mass,...
        'alf',sfrAlfa,'afac',sfrAfac);
    sfrProf(:,i)=sfr;
    dmass=sfr.*(dt.*1e9);
    
    stellarMass(i)=stellarMass(i-1)+sum(dmass);
    
    % update the gas mass.
    gal.gasMass.mass=max(gal.gasMass.mass-dmass,zeros(size(dmass)));
    gasProf(:,i)=gal.gasMass.mass;
    
    gasMass(i)=sum(gasProf(:,i));
    
end

if wbarFlag; close(hb) ; end


res.time=time;
res.rpos=rp;
res.pram=pram0;
res.rhoICM=rhoICM;

res.rgal=gal.gasMass.rr(2:end);
res.gasProf=gasProf;
res.sfrProf=sfrProf;
res.sfr=sum(sfrProf,1);
res.ssfr=res.sfr./stellarMass;
res.gasMass=sum(gasProf,1);
res.stellarMass=stellarMass;


if plotFlag
    
    
    
    %% figures
    
    timeLab='Time [Gyr]';
    % plot mstar & mgas
    
    cmap=brewermap(8,'Set1');
    
    h=[];
    figure
    h(1)=plot(time,mgas(1:ind)./mgas(1),'color',cmap(1,:),'DisplayName','SF Only');
    hold on
    h(2)=plot(time,gasMass./gasMass(1),'color',cmap(2,:),'DisplayName','SF \& RPS');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
    
    
    grid
    xlabelmine(timLab);
    ylabelmine('Gas Mass Loss');%$[\mathrm{M_\odot}]$')
    set(gca,'Fontsize',14)
    
    %printout_fig(gcf,['sfRps_' orbTag '_compare_mgas'],'v')
    
    % ------------------------------
    
    h=[];
    figure
    h(1)=plot(time,mstar(1:ind)./mstar(1),'color',cmap(1,:),'DisplayName','SF Only');
    hold on
    h(2)=plot(time,stellarMass./stellarMass(1),'color',cmap(2,:),'DisplayName','SF \& RPS');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14,'location','SouthEast')
    
    
    grid
    xlabelmine(timLab);
    ylabelmine('Stellar Mass Growth'); %$[\mathrm{M_\odot}]$')
    set(gca,'Fontsize',14)
    
    %printout_fig(gcf,['sfRps_' orbTag '_compare_mstar'],'v')
    
    % ------------------------------
    
    h=[];
    figure
    h(1)=plot(time,sum(sfr1(:,1:ind),1)./sum(sfr1(:,1),1),'color',cmap(1,:),'DisplayName','SF Only');
    hold on
    h(2)=plot(time,sum(sfrProf,1)./sum(sfrProf(:,1),1),'color',cmap(2,:),'DisplayName','SF \& RPS');
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
    
    
    grid
    xlabelmine(timeLab);
    ylabelmine('SFR Decay'); %$[\mathrm{M_\odot\,yr^{-1}}]$')
    set(gca,'Fontsize',14)
    
    %printout_fig(gcf,['sfRps_' orbTag '_compare_sfr'],'v')
    
    % ------------------------------
    
    h=[];
    figure
    h(1)=plot(time,sum(sfr1(:,1:ind),1)./mstar(1:ind),'color',cmap(1,:),'DisplayName','SF Only');
    hold on
    h(2)=plot(time,sum(sfrProf,1)./stellarMass,'color',cmap(2,:),'DisplayName','SF \& RPS');
    plot(time,ones(size(time)).*1e-11,'--k')
    
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
    
    
    grid
    xlabelmine(timeLab);
    ylabelmine('sSFR $[\mathrm{yr^{-1}}]$')
    set(gca,'Fontsize',14)
    
    %printout_fig(gcf,['sfRps_' orbTag '_compare_ssfr'],'v')
    
    % ------------------------------
    
    db=floor(ind/12);
    ii=1:db:ind;
    ii(end)=ind;
    
    figure
    h=[];
    cnt=0;
    for i=ii
        cnt=cnt+1;
        tag=sprintf('%3.2f Gyr',time(i));
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
    
    %printout_fig(gcf,['sf_' orbTag '_gasProf'],'v')
    
    figure
    h=[];
    cnt=0;
    for i=ii
        cnt=cnt+1;
        tag=sprintf('%3.2f Gyr',time(i));
        h(cnt)=plot(gal.gasMass.rr(2:end),gasProf(:,i),'Displayname',tag);
        hold on
    end
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
    
    grid
    xlabelmine('galactic radius [kpc]');
    ylabelmine('Gas Mass $[\mathrm{M_\odot}]$')
    set(gca,'Fontsize',14)
    titlemine('SF \& RPS')
    
    %printout_fig(gcf,['sfRPS_' orbTag '_gasProf'],'v')
    
    
    % -----------------------
    
    
    figure
    h=[];
    cnt=0;
    for i=ii
        cnt=cnt+1;
        tag=sprintf('%3.2f Gyr',time(i));
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
    
    %printout_fig(gcf,['sf_' orbTag '_sfrProf'],'v')
    
    figure
    h=[];
    cnt=0;
    for i=ii
        cnt=cnt+1;
        tag=sprintf('%3.2f Gyr',time(i));
        h(cnt)=plot(gal.gasMass.rr(2:end),sfrProf(:,i),'Displayname',tag);
        hold on
    end
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
    
    grid
    xlabelmine('galactic radius [kpc]');
    ylabelmine('SFR $[\mathrm{M_\odot \,yr^{-1}}]$')
    set(gca,'Fontsize',14)
    titlemine('SF \& RPS')
    
    %%printout_fig(gcf,['sfRPS_' orbTag '_sfrProf'],'v')
    
    
    
end




