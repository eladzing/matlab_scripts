function [res, resBird]=get_halo_history_profiles(galID,varargin)
%% we want to follow a galaxy back over time and record some important information.

% default
global illUnits
global cosmoStruct
%global fullSnapMask
global BASEPATH

zmax=2.;
snapStart=99;
verboseFlag=false;
treeType='SubLink_gal';
nbins=400;
%    units; % load general unit structure in cgs.

ii=1;
while(ii<=length(varargin))
    switch(lower(varargin{ii}))
        case {'zmax','maxz'}
            ii=ii+1;
            zmax=varargin{ii};
        case {'snapstart','startsnap'}
            ii=ii+1;
            snapStart=varargin{ii};
        case {'treetype','type','tree'}
            ii=ii+1;
            treeType=varargin{ii};
        case {'verbose','v'}
            verboseFlag=true;
        case {'nbins','nb'}
            ii=ii+1;
            nbins=varargin{ii};
        otherwise
            error('GET_HALO_HISTORY_PROFILES - Illegal argument: %s',varargin{ii});
    end
    ii=ii+1;
end

% get the tree

treeFields={'SnapNum','SubfindID','SubhaloVel','SubhaloPos','SubhaloHalfmassRadType','SubhaloGrNr',...
    'SubhaloMassInRadType','SubhaloSFRinRad','Group_M_Crit200','SubfindID','SubhaloHalfmassRadType',...
    'Group_R_Crit200','GroupPos'};
tree = illustris.sublink.loadTree(BASEPATH,snapStart,galID,treeFields,true,treeType);

% start climbing the tree

redshifts=illustris.utils.snap2redshift(tree.SnapNum);
[~,lastInd]=min(abs(redshifts-zmax));


% get unit conversion factors 
unitFactors=illustris.utils.calc_illUnits(tree.SnapNum(1:lastInd));

%firstInd=length(tree.SnapNum);

%%
res.Mvir=tree.Group_M_Crit200(1:lastInd)'.*unitFactors.massUnit;
res.Rvir=tree.Group_R_Crit200(1:lastInd)'.*unitFactors.lengthUnit;
[~, ~, tv, ~]=calculate_virials('mvir',res.Mvir,'cosmo',cosmoStruct,'delv',200,'crit');
res.Tvir=tv; 

res.starHalfMassRad=tree.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,1:lastInd).*unitFactors.lengthUnit;
res.gasalfMassRad=tree.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,1:lastInd).*unitFactors.lengthUnit;

res.SubfindID=tree.SubfindID(1:lastInd)';
res.sfr=tree.SubhaloSFRinRad(1:lastInd)';
res.stellarMass=tree.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,1:lastInd).*unitFactors.massUnit;
res.ssfr=res.sfr./res.stellarMass+...
    10^0.5*1e-17.*10.^(0.5.*rand(size(res.sfr)));

%aexp=illustris.utils.snap2aexpn(tree.SnapNum(1:lastInd));
res.zred=illustris.utils.snap2redshift(tree.SnapNum(1:lastInd));
%res.zredExt=res.zred(fullSnapMask(tree.SnapNum(1:lastInd)+1));

%cntExt=0;


%% set to zero
%len=lastInd;
fofID=tree.SubhaloGrNr(1:lastInd);
%% climb up the tree

for i=1:lastInd
    if verboseFlag
        perc=floor(100*i/lastInd);
        if mod(perc,5)==0
            fprintf('completed %i %% of tree climbing  \n',floor(100*i/lastInd));
        end
    end
    % read information from the the relevant snapshot
    snap=tree.SnapNum(i);
    
    illustris.utils.set_illUnits(snap) % set the simulation units for the snapshot.
    
    
    %     % load gas from in sub halo
    %     if fullSnapMask(snap+1)
    %         cntExt=cntExt+1;
    %         gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','Velocities','StarFormationRate',...
    %             'GFM_CoolingRate','EnergyDissipation','Machnumber'};
    %     else
    %           end
    
    gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','Velocities','StarFormationRate'};
    gas=illustris.snapshot.loadSubhalo(BASEPATH, tree.SnapNum(i), fofID(i), 'gas',gasFields);
    
    if gas.count==0
        continue
    end
    
    %sfMask=gas.StarFormationRate==0;
    
    nDensity=gas.Density.*illUnits.numberDensityFactor; %in cm^-3
    gas=illustris.utils.addTemperature( gas );
    gas=illustris.utils.addEntropy( gas );
    gas=illustris.utils.addPressure( gas );
    mass=gas.Masses.*illUnits.massUnit; %in Solarmass
    
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,tree.GroupPos(:,i)).*illUnits.lengthUnit;
    gas=illustris.utils.addCellRadius(gas);
    
    mask=gas.StarFormationRate==0;
    %gas=illustris.utils.addCellRadius(gas);
    
    
    
    
    
    % deal with velocitites
%     for k=1:3
%         vv=gas.Velocities(k,:).*illustris.utils.velocityFactor(tree.SnapNum(i),'gas')-...
%             tree.SubhaloVel(k,i).*illustris.utils.velocityFactor(tree.SnapNum(i),'sub');
%     end
%     vel=sqrt(sum(vv.^2,1));
%     clear vv
    
    
    
    % find distance from galaxy center in comoving coordinates !
%    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,tree.SubhaloPos(:,i));
%    gasDist=sqrt( sum(gas.newCoord.^2,1));
    
%     rhalfStar=tree.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,i); % stellar half mass radius
%     rhalfGas=tree.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,i); % gas half mass radius
    
    
%     % generate masks for diffrent gaseous phase components.
%     polys=phaseDiagram_polygons;
%     
%     
%     %sfMask          =inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
%     sfMask          =gas.StarFormationRate>0;
%     coldDenseMask   = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.coldDense_polygonMask(:,1),polys.coldDense_polygonMask(:,2));
%     coldDiluteDMask = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.coldDilute_polygonMask(:,1),polys.coldDilute_polygonMask(:,2));
%     warmHotMask     = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.warmHot_polygonMask(:,1),polys.warmHot_polygonMask(:,2));
%     hotMask         = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.hotTemp_polygonMask(:,1),polys.hotTemp_polygonMask(:,2));
%     
%     
%     if fullSnapMask(snap+1)
%         tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr^-1
%     end
%     
%     
    %% calculate profiles 
      % Temperature 
    res.profT(i) = mk_radial_profile_cells(gas.newCoord,gas.CellRadius,gas.Temperature,'wt',mass,'type','intensive','nb',nbins,'mask',mask);
    
    % Entropy
    res.profK(i) = mk_radial_profile_cells(gas.newCoord,gas.CellRadius,gas.Entropy,'wt',mass,'type','intensive','nb',nbins,'mask',mask);
    % Pressure
    res.profP(i) = mk_radial_profile_cells(gas.newCoord,gas.CellRadius,gas.Pressure,'wt',mass,'type','intensive','nb',nbins,'mask',mask);
%     res.profP2(i) = mk_radial_profile_cells(gas.newCoord,gas.CellRadius,gas.Pressure,'wt',ones(size(mass)),'type','intensive','nb',nbins,'mask',mask);
%     res.profP3(i) = mk_radial_profile_cells(gas.newCoord,gas.CellRadius,gas.Pressure,'wt',gas.CellRadius,'type','intensive','nb',nbins,'mask',mask);
    res.profM(i)  = mk_radial_profile_cells(gas.newCoord,gas.CellRadius,mass,'wt',ones(size(mass)),'type','extensive','nb',nbins,'mask',mask);
    
    
    
    % Density
end

resBird=0;
    
end
    
    
    
    
    




