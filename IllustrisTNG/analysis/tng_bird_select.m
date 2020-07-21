

%% create parameter plot which defines selection

%% read in gas
if readFlag
    id=0;
    gas=illustris.snapshot.loadSubhalo(bp, 99, id, 'gas',...
        {'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses'});
    
end
if ~exist('firstPassFlag','var')
    firstPassFlag=true;
end

%% set limits
if firstPassFlag
    %id=1; % index in subs structure (0-based)
    
    tlim=[3.9 8.1];
    nlim=[-6 1];
    
    rhoFac=densityUnit.*(Units.Ms/Units.kpc^3/Units.mm); %in cm^-3
    
    dens=gas.Density.*rhoFac; %in cm^-3
    if ~isfield(gas,'Temperature')
        gas.Temperature=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
    end
    temp=gas.Temperature; %   illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
    mass=gas.Masses.*massUnit; %in Solarmass
    cellSizes=(gas.Masses./gas.Density).^(1/3);
    
    
    
    rhalfStar= subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1); % stellar half mass radius
    
    %% define coordinates, and boxsize
    
    if ~isfield(gas,'newCoord')
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
    end
    
    dim=find(size(gas.newCoord)~=3);
    m1=abs(max(gas.newCoord,[],dim));
    m2=abs(min(gas.newCoord,[],dim));
    ll=max(cat(1,m1,m2));
    
    boxx=ceil(2.*ll);
    clear m1 m2 ll
    
    %%  distances and mask
    gasDist=sqrt( sum(gas.newCoord.^2,1));
    distMask=true(size(mass));%    gasDist<=2.0.*rhalf;
    
    
    %% get birds
    [bird, binsize, xxlim,yylim]= histogram2d(log10(temp(distMask)),log10(dens(distMask)),mass(distMask),...
        'xxlim',tlim,'yylim',nlim,'len',256-1);
    
    [birdRT, ~, xxlimRT,yylimRT]= histogram2d(gasDist(distMask),log10(temp(distMask)),mass(distMask),...
        'yylim',tlim,'len',256-1);
    
    [birdNT, ~, xxlimNT,yylimRN]= histogram2d(gasDist(distMask),log10(dens(distMask)),mass(distMask),...
        'yylim',nlim,'len',256-1);
    
    
    
end

%% plot birds
hf=figure('pos',[100 100 2000 1200]);
ah=subplot(2,2,1);
birdP=bird(:,:,1)./sum(sum(bird(:,:,1)));
illustris.plots.plot_phaseDiagram(xxlim,yylim,birdP,'lines','clims',[-7 0],'axes',ah,'figure',hf,'nt','title','')


cmap=brewermap(256,'*Spectral');
cmap(1,:)=[1 1 1];
colormap(cmap);

subplot(2,2,2);
imagesc(xxlimRT,yylimRT,log10(squeeze(birdRT(:,:,1))./sum(birdRT(:))))
set(gca,'Ydir','normal','Fontsize',12)
xlabelmine('$r\,[\mathrm{kpc/h}]$')
ylabelmine('$\log(T)\,[\mathrm{K}]$')
hold on
plot(2*rhalfStar.*[1 1],yylimRT,'--k');
grid

subplot(2,2,4);
imagesc(xxlimNT,yylimRN,log10(squeeze(birdNT(:,:,1))./sum(birdNT(:))))
set(gca,'Ydir','normal','Fontsize',12)
xlabelmine('$r\,[\mathrm{kpc/h}]$')
ylabelmine('$\log(n)\,[\mathrm{cm^{-3}}]$')
hold on
plot(2*rhalfStar.*[1 1],yylimRN,'--k');
grid




%% select points to map
[xv, yv] = getline(ah); % uncomment this line for interactive selection
line(xv,yv, 'Color','k');

dens_log = log10(dens);
temp_log = log10(temp);

birdMask = inpolygon(dens_log,temp_log,xv,yv);


totMask=distMask & birdMask;

%% creat t & ro vs distance

[birdRT2, binsizeRT, xxlimRT2,yylimRT2]= histogram2d(gasDist(totMask),log10(temp(totMask)),mass(totMask),...
    'yylim',tlim,'len',256-1);

[birdNT2, binsizeNT, xxlimNT2,yylimNT2]= histogram2d(gasDist(totMask),log10(dens(totMask)),mass(totMask),...
    'yylim',nlim,'len',256-1);


%% plot the maps

% create uniform grid
%if firstPassFlag
ng=256;
cub=cell2grid(gas.newCoord(:,totMask),gas.Masses(totMask),cellSizes(totMask),'extensive','boxSide',boxx,'Ngrid',ng);
%end

circ.radius=2.*rhalfStar;
circ.color='k';


hf2=figure('pos',[100 100 2000 1200]);
ah=subplot(2,4,1);

% plot zy
subplot(4,2,[6 8])
illustris.plots.mkmapGas('data',cub,'yz','normfactor',rhoFac,'log','thick',-ng,'axes',gca,'fig',gcf,'grid','circ',circ)


% plot xy
subplot(4,2,[5 7])
illustris.plots.mkmapGas('data',cub,'xy','normfactor',rhoFac,'log','thick',-ng,'axes',gca,'fig',gcf,'grid','circ',circ)

% plot zx
subplot(4,2,[1 3])
illustris.plots.mkmapGas('data',cub,'zx','normfactor',rhoFac,'log','thick',-ng,'axes',gca,'fig',gcf,'grid','circ',circ)

subplot(4,2,2)
imagesc(xxlimRT2,yylimRT2,log10(squeeze(birdRT2(:,:,1))./sum(birdRT2(:))))
set(gca,'Ydir','normal','Fontsize',12)
xlabelmine('$r\,[\mathrm{kpc/h}]$')
ylabelmine('$\log(T)\,[\mathrm{K}]$')
hold on
plot(2*rhalfStar.*[1 1],yylimRN,'--k');
grid

subplot(4,2,4);
imagesc(xxlimNT2,yylimNT2,log10(squeeze(birdNT2(:,:,1))./sum(birdNT2(:))))
set(gca,'Ydir','normal','Fontsize',12)
xlabelmine('$r\,[\mathrm{kpc/h}]$')
ylabelmine('$\log(n)\,[\mathrm{cm^{-3}}]$')
hold on
plot(2*rhalfStar.*[1 1],yylimRN,'--k');
grid



