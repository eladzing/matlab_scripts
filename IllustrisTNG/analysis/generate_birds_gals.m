%% Study the differences between the gas in sf vs. quiescent centals
% go over all centrals and extract phase-space plots for the gas in
% several areas: within 2*r_h, all gas in the subhalo, others?

%% set framework

%global matFilePath
global cosmoStruct
global massUnit
global densityUnit
%global hub
sim='100';
snap=99; %z=0
bp=illustris.set_env(sim,'draco');

if ~exist('readFlag')
    readFlag=true;
end

if readFlag
    fprintf(' *** Reading data *** \n');
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
end



massThresh=10^9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
units; % load general unit structure in cgs.
%% selec galaxies

massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*massUnit; % stellar mass within 2*rhalf
massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
mask= subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas;

% create index list of central galaxies
galList=0:subs.count-1;
galList=galList(mask);


%% get polygon masks 
plotPolygonFlag=false;

phaseDiagram_polygons   %% call script defining polygons. 



%% generate birds
fprintf(' *** Running over %s Galaxies *** \n',num2str(length(galList)))

%% set limits
tlim=[3.9 8.1];
nlim=[-6 1];

cnt=0;
step=5;
stepNext=5;
for id=galList
    cnt=cnt+1;
    galBirds(cnt).id=id;
    
    
    perCent=floor(cnt/length(galList).*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    % load gas from in sub halo
    gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas', ...
        {'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses'});
    
    % get density and temperature
    dens=gas.Density.*densityUnit.*(Units.Ms/Units.kpc^3/Units.mm); %in cm^-3
    temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
    mass=gas.Masses.*massUnit; %in Solarmass
    
    
    % get virial temperatures 
    [~, ~, tvirHost, ~]=calculate_virials('mvir',fofs.Group_M_Mean200(subsInfo.hostFof(id+1)+1)*massUnit,'cosmo',cosmoStruct,'delv',200);
    [~, ~, tvirSubs, ~]=calculate_virials('mvir',subs.SubhaloMass(id+1)*massUnit,'cosmo',cosmoStruct,'delv',200);

    
    % find distance from galaxy center
    
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
    gasDist=sqrt( sum(gas.newCoord.^2,1));
    %gasDist=sqrt( (gas.Coordinates(1,:)-subs.SubhaloPos(1,id)).^2+...
    %    (gas.Coordinates(2,:)-subs.SubhaloPos(2,id)).^2+...
    %    (gas.Coordinates(3,:)-subs.SubhaloPos(3,id)).^2);
    
    
    
    % find importat radii
    rmax=max(gasDist);
    rhalfStar=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1); % stellar half mass radius
    rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1); % stellar half mass radius
    galBirds(cnt).rmax=rmax;
    galBirds(cnt).rhalfStar=rhalfStar;
    galBirds(cnt).rhalfGas=rhalfGas;
    
    
       
    %% for gas within the gal (2*r_half - stellar)
    
    mask=gasDist<=2.0.*rhalfStar;
    
   
    
    if any(mask)
        
        ro=dens(mask);
        tt=temp(mask);
        mm=mass(mask);
        
        % [bird, binsize, xxlim,yylim]= histogram2d(log10(mass(mask)),gr300(mask),ones(size(gr300(mask))),'xxlim',[9 12],'yylim',[0 1],'len',200);
        [bird1, binsize1, xxlim1,yylim1]= histogram2d(log10(tt),log10(ro),mm,...
            'xxlim',tlim,'yylim',nlim,'len',200-1);
        
        galBirds(cnt).BirdInGal=bird1;
        
        galBirds(cnt).GasInGal.total=sum(mm);
        galBirds(cnt).GasInGal.starForming=sum(mm(inpolygon(log10(ro),log10(tt),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2))));
        galBirds(cnt).GasInGal.coldDense=sum(mm(inpolygon(log10(ro),log10(tt),polys.coldDense_polygonMask(:,1),polys.coldDense_polygonMask(:,2))));
    
        galBirds(cnt).GasInGal.coldDilute=sum(mm(inpolygon(log10(ro),log10(tt),polys.coldDilute_polygonMask(:,1),polys.coldDilute_polygonMask(:,2))));
        galBirds(cnt).GasInGal.warmHot=sum(mm(inpolygon(log10(ro),log10(tt),polys.warmHot_polygonMask(:,1),polys.warmHot_polygonMask(:,2))));
        galBirds(cnt).GasInGal.hotTemp=sum(mm(inpolygon(log10(ro),log10(tt),polys.hotTemp_polygonMask(:,1),polys.hotTemp_polygonMask(:,2))));
         
        galBirds(cnt).GasInGal.aboveTvirSub=sum(mm(tt>tvirSubs   ));
        galBirds(cnt).GasInGal.aboveTvirHost=sum(mm(tt>tvirHost ));
    
        
    else
        galBirds(cnt).BirdInGal=0;
        galBirds(cnt).GasInGal.total=0;
        galBirds(cnt).GasInGal.starForming=0;
        galBirds(cnt).GasInGal.coldDense=0;
    
        galBirds(cnt).GasInGal.coldDilute=0;
        galBirds(cnt).GasInGal.warmHot=0;
        galBirds(cnt).GasInGal.hotTemp=0;
         
        galBirds(cnt).GasInGal.aboveTvirSub=0;
        galBirds(cnt).GasInGal.aboveTvirHost=0;
    end
    
    
   
    
    
    
    %% for gas within the cgm (2*r_half  - r_half gas)
    
    mask=gasDist>2.0.*rhalfStar &  ...
        gasDist<=rhalfGas;
    
   
    
    if any(mask)
        
         ro=dens(mask);
    tt=temp(mask);
    mm=mass(mask);
    
    
    [bird2, binsize2, xxlim2,yylim2]= histogram2d(log10(tt),log10(ro),mm,...
        'xxlim',tlim,'yylim',nlim,'len',200-1);
    tr
        galBirds(cnt).birdInCGM=bird2;
    
        
        galBirds(cnt).GasInCGM.total=sum(mm);
        galBirds(cnt).GasInCGM.starForming=sum(mm(inpolygon(log10(ro),log10(tt),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2))));
        galBirds(cnt).GasInCGM.coldDense=sum(mm(inpolygon(log10(ro),log10(tt),polys.coldDense_polygonMask(:,1),polys.coldDense_polygonMask(:,2))));
    
        galBirds(cnt).GasInCGM.coldDilute=sum(mm(inpolygon(log10(ro),log10(tt),polys.coldDilute_polygonMask(:,1),polys.coldDilute_polygonMask(:,2))));
        galBirds(cnt).GasInCGM.warmHot=sum(mm(inpolygon(log10(ro),log10(tt),polys.warmHot_polygonMask(:,1),polys.warmHot_polygonMask(:,2))));
        galBirds(cnt).GasInCGM.hotTemp=sum(mm(inpolygon(log10(ro),log10(tt),polys.hotTemp_polygonMask(:,1),polys.hotTemp_polygonMask(:,2))));
         
        galBirds(cnt).GasInCGM.aboveTvirSub=sum(mm(tt>tvirSubs    ));
        galBirds(cnt).GasInCGM.aboveTvirHost=sum(mm(tt>tvirHost ));
    
        
    else
        galBirds(cnt).BirdInCGM=0;
        galBirds(cnt).GasInCGM.total=0;
        galBirds(cnt).GasInCGM.starForming=0;
        galBirds(cnt).GasInCGM.coldDense=0;
    
        galBirds(cnt).GasInCGM.coldDilute=0;
        galBirds(cnt).GasInCGM.warmHot=0;
        galBirds(cnt).GasInCGM.hotTemp=0;
         
        galBirds(cnt).GasInCGM.aboveTvirSub=0;
        galBirds(cnt).GasInCGM.aboveTvirHost=0;
           
    end
    
    
    %% gas within 2* half stellar mass radius and edge
    mask=gasDist>2.0.*rhalfStar ; 
    
    if any(mask)
        
         ro=dens(mask);
    tt=temp(mask);
    mm=mass(mask);
    
    
    [bird3, binsize3, xxlim3,yylim3]= histogram2d(log10(tt),log10(ro),mm,...
        'xxlim',tlim,'yylim',nlim,'len',200-1);
    
    galBirds(cnt).birdInSub=bird3;
    
    galBirds(cnt).GasInSub.total=sum(mm);
        galBirds(cnt).GasInSub.starForming=sum(mm(inpolygon(log10(ro),log10(tt),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2))));
        galBirds(cnt).GasInSub.coldDense=sum(mm(inpolygon(log10(ro),log10(tt),polys.coldDense_polygonMask(:,1),polys.coldDense_polygonMask(:,2))));
    
        galBirds(cnt).GasInSub.coldDilute=sum(mm(inpolygon(log10(ro),log10(tt),polys.coldDilute_polygonMask(:,1),polys.coldDilute_polygonMask(:,2))));
        galBirds(cnt).GasInSub.warmHot=sum(mm(inpolygon(log10(ro),log10(tt),polys.warmHot_polygonMask(:,1),polys.warmHot_polygonMask(:,2))));
        galBirds(cnt).GasInSub.hotTemp=sum(mm(inpolygon(log10(ro),log10(tt),polys.hotTemp_polygonMask(:,1),polys.hotTemp_polygonMask(:,2))));
         
        galBirds(cnt).GasInSub.aboveTvirSub=sum(mm(tt>tvirSubs));
        galBirds(cnt).GasInSub.aboveTvirHost=sum(mm(tt>tvirHost));
    
        
    else
        galBirds(cnt).BirdInSub=0;
        galBirds(cnt).GasInSub.total=0;
        galBirds(cnt).GasInSub.starForming=0;
        galBirds(cnt).GasInSub.coldDense=0;
    
        galBirds(cnt).GasInSub.coldDilute=0;
        galBirds(cnt).GasInSub.warmHot=0;
        galBirds(cnt).GasInSub.hotTemp=0;
         
        galBirds(cnt).GasInSub.aboveTvirSub=0;
        galBirds(cnt).GasInSub.aboveTvirHost=0;
           
    end
        

    
    
    %% add additional info to struct
    
    galBirds(cnt).stellarMass=massAllGals(id+1);
    galBirds(cnt).xxlim=xxlim3;
    galBirds(cnt).yylim=yylim3;
    galBirds(cnt).binsize=binsize3;
    galBirds(cnt).TvirSub=tvirSubs;
    galBirds(cnt).TvirHost=tvirHost;
    
    
 
    
    
    
    
end


%% save to mat file
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    fname=sprintf('galaxy_birds_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'galBirds','cnt','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
end
fprintf(' *** DONE!  *** \n');


