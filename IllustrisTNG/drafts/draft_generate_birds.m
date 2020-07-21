%% Study the differences between the gas in sf vs. quiescent centals
% go over all centrals and extract phase-space plots for the gas in
% several areas: within 2*r_h, all gas in the subhalo, others?

%% set framework

%global matFilePath
global massUnit
global densityUnit
global hub
sim='100';
snap=99; %z=0
% bp=illustris.set_env(sim);
% 
% fofs=illustris.groupcat.loadHalos(bp,snap);
% subs=illustris.groupcat.loadSubhalos(bp,snap);
% 


massThresh=1e9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
units; % load general unit structure in cgs.
%% selec galaxies

massAllGals= subs.SubhaloMassInRadType(5,:).*massUnit; % stellar mass within 2*rhalf
massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
mask=subsInfo.hasDM & subsInfo.hasStars & massMask & subsInfo.hostHasVirial;

% create index list of central galaxies
galList=1:subs.count;
galList=galList(mask);


%% generate birds
cnt=0;
for id=galList
    cnt=cnt+1;    
    galBirds(cnt).id=id;
    

    % load gas from in sub halo
    gas=illustris.snapshot.loadSubhalo(bp, snap, id, illustris.partTypeNum('gas'), ...
        {'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses'});
    
    % get density and temperature
    dens=gas.Density.*densityUnit.*(Units.Ms/Units.kpc^3/Units.mm); %in cm^-3
    temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
    mass=gas.Masses.*massUnit; %in Solarmass
    
    % find distance from galaxy center
    
    newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id));
    gasDist=sqrt( sum(newCoord.^2,1));
    %gasDist=sqrt( (gas.Coordinates(1,:)-subs.SubhaloPos(1,id)).^2+...
    %    (gas.Coordinates(2,:)-subs.SubhaloPos(2,id)).^2+...
    %    (gas.Coordinates(3,:)-subs.SubhaloPos(3,id)).^2);
    
    
    
    % find importnat radii
    rmax=max(gasDist);
    rhalf=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id); % stellar half mass radius
    galBirds(cnt).rmax=rmax;
    galBirds(cnt).rhalf=rhalf;
    
    %% set limits 
    tlim=[3.9 8.1];
    nlim=[-4.1 0.1];
    
    %% for gas within the gal (2*r_half - stellar)
    
    mask=gasDist<=2.0.*rhalf;
    
    % [bird, binsize, xxlim,yylim]= histogram2d(log10(mass(mask)),gr300(mask),ones(size(gr300(mask))),'xxlim',[9 12],'yylim',[0 1],'len',200);
    [bird1, binsize1, xxlim1,yylim1]= histogram2d(log10(temp(mask)),log10(dens(mask)),mass(mask),...
        'xxlim',tlim,'yylim',nlim,'len',200-1);
    
    galBirds(cnt).GasInGal=sum(mass(mask));
    galBirds(cnt).BirdInGal=bird1;
    
    
    %% for gas within the cgm (2*r_half  - 0.5
    
    mask=gasDist>2.0.*rhalf &  ...
        gasDist<=0.5.*rmax;
    
    [bird2, binsize2, xxlim2,yylim2]= histogram2d(log10(temp(mask)),log10(dens(mask)),mass(mask),...
            'xxlim',tlim,'yylim',nlim,'len',200-1);
    
    galBirds(cnt).GasInCGM=sum(mass(mask));
    galBirds(cnt).birdInCGM=bird2;
    
    %% gas within 0.5rmax and rmax
    mask=gasDist>2.0.*rhalf & gasDist <=rmax;
    
    [bird3, binsize3, xxlim3,yylim3]= histogram2d(log10(temp(mask)),log10(dens(mask)),mass(mask),...
            'xxlim',tlim,'yylim',nlim,'len',200-1);
    
    galBirds(cnt).GasInSub=sum(mass(mask));
    galBirds(cnt).birdInSub=bird3;
    
    
    galBirds(cnt).stellarMass=massAllGals(id);
    
    
end