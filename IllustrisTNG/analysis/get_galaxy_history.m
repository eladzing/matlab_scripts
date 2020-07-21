function res=get_galaxy_history(galID,varargin)
%% we want to follow a galaxy back over time and record some important information.

% default


global illUnits
global fullSnapMask
global BASEPATH

zmax=2.;
snapStart=99;
verboseFlag=false;
treeType='SubLink_gal';
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
        otherwise
            error('GET_GALAXY_HISTORY - Illegal argument: %s',varargin{ii});
    end
    ii=ii+1;
end



% get the tree
treeFields={'SnapNum','SubfindID','SubhaloVel','SubhaloPos','SubhaloHalfmassRadType',...
    'SubhaloMassInRadType','SubhaloSFRinRad'};
tree = illustris.sublink.loadTree(BASEPATH,snapStart,galID,treeFields,true,treeType);




% start climbing the tree

redshifts=illustris.utils.snap2redshift(tree.SnapNum);
[~,lastInd]=min(abs(redshifts-zmax));
%firstInd=length(tree.SnapNum);

%%
res.sfr=tree.SubhaloSFRinRad(1:lastInd)';
res.stellarMass=tree.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,1:lastInd).*illUnits.massUnit;
res.ssfr=res.sfr./res.stellarMass+...
    10^0.5*1e-17.*10.^(0.5.*rand(size(res.sfr)));

aexp=illustris.utils.snap2aexpn(tree.SnapNum(1:lastInd));
res.zred=illustris.utils.snap2redshift(tree.SnapNum(1:lastInd));
res.zredExt=res.zred(fullSnapMask(tree.SnapNum(1:lastInd)+1));

cntExt=0;


%% set to zero
len=lastInd;
lenExt=length(res.zredExt);

inGal.meanTempMW=zeros(2,len);
inGal.medianTemp=zeros(2,len);
inGal.meanEntMW=zeros(2,len);
inGal.medianEnt=zeros(2,len);
inGal.meanDensN=zeros(2,len);
inGal.medianDensN=zeros(2,len);
inGal.velDisp=zeros(2,len);
inGal.velDispMW=zeros(2,len);
inGal.gasMass=zeros(2,len);
inGal.avgDensN=zeros(2,len);
inGal.cellNum=zeros(2,len);

% data found only in full snapshots.
inGal.meanTcMW=zeros(2,lenExt);
inGal.medianTc=zeros(2,lenExt);
inGal.EnergyDissipation=zeros(2,lenExt);
inGal.meanMachEW=zeros(2,lenExt);
inGal.meanMach=zeros(2,lenExt);
inGal.medianMach=zeros(2,lenExt);


% mass in different gas components

% star forming gas
inGal.sfGas.mass=zeros(1,len);
inGal.sfGas.EnergyDissipation=zeros(1,lenExt);
inGal.sfGas.meanMachEW=zeros(1,lenExt);

inGal.coldDenseGas=inGal.sfGas;
inGal.coldDiluteGas=inGal.sfGas;
inGal.warmHotGas=inGal.sfGas;
inGal.hotGas=inGal.sfGas;

inCGM=inGal;
inSub=inGal;


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
    
    
    % load gas from in sub halo
    if fullSnapMask(snap+1)
        cntExt=cntExt+1;
        gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','Velocities','StarFormationRate',...
            'GFM_CoolingRate','EnergyDissipation','Machnumber'};
    else
        gasFields={'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','Velocities','StarFormationRate'};
    end
    
    
    gas=illustris.snapshot.loadSubhalo(BASEPATH, tree.SnapNum(i), tree.SubfindID(i), 'gas',gasFields);
    
    if gas.count==0
        continue
    end
    
    
    nDensity=gas.Density.*illUnits.numberDensityFactor; %in cm^-3
    gas=illustris.utils.addTemperature( gas );
    gas=illustris.utils.addEntropy( gas );
    mass=gas.Masses.*illUnits.massUnit; %in Solarmass
    
    % deal with velocitites
    for k=1:3
        vv=gas.Velocities(k,:).*illustris.utils.velocityFactor(tree.SnapNum(i),'gas')-...
            tree.SubhaloVel(k,i).*illustris.utils.velocityFactor(tree.SnapNum(i),'sub');
    end
    vel=sqrt(sum(vv.^2,1));
    clear vv
    
    
    
    % find distance from galaxy center in comoving coordinates !
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,tree.SubhaloPos(:,i));
    gasDist=sqrt( sum(gas.newCoord.^2,1));
    
    rhalfStar=tree.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,i); % stellar half mass radius
    rhalfGas=tree.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,i); % gas half mass radius
    
    
    % generate masks for diffrent gaseous phase components.
    polys=phaseDiagram_polygons;
    
    
    %sfMask          =inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
    sfMask          =gas.StarFormationRate>0;
    coldDenseMask   = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.coldDense_polygonMask(:,1),polys.coldDense_polygonMask(:,2));
    coldDiluteDMask = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.coldDilute_polygonMask(:,1),polys.coldDilute_polygonMask(:,2));
    warmHotMask     = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.warmHot_polygonMask(:,1),polys.warmHot_polygonMask(:,2));
    hotMask         = ~sfMask & inpolygon(log10(nDensity),log10(gas.Temperature),polys.hotTemp_polygonMask(:,1),polys.hotTemp_polygonMask(:,2));
    
    
    if fullSnapMask(snap+1)
        tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr
    end
    
    
    % generate the relevant information
    
    
    
    %% for gas within the gal (2*r_half - stellar)
    
    distMask=gasDist<=2.0.*rhalfStar;
    
    for k=1:2
        switch k
            case 1
                mask=distMask & ~sfMask;
            case 2
                mask=distMask;
        end
        
        
        
        if any(mask)
            mm=mass(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            vv=vel(mask);
            nDens=nDensity(mask);
            
            % calc bird
            
            inGal.meanTempMW(k,i)=sum(mm.*tmp)/sum(mm);
            inGal.medianTemp(k,i)=median(tmp);
            %                 inGal.quant25Temp(k,i)=quantile(tmp,0.25);
            %                 inGal.quant75Temp(k,i)=quantile(tmp,0.75);
            %                 inGal.stdTemp(k,i)=std(tmp);
            
            inGal.meanEntMW(k,i)=sum(mm.*ent)/sum(mm);
            inGal.medianEnt(k,i)=median(ent);
            %                 inGal.quant25Ent(k,i)=quantile(ent,0.25);
            %                 inGal.quant75Ent(k,i)=quantile(ent,0.75);
            %                 inGal.stdEnt(k,i)=std(ent);
            
            inGal.meanDensN(k,i)=mean(nDens);
            inGal.medianDensN(k,i)=median(nDens);
            %                 inGal.quant25DensN=quantile(nDens,0.25);
            %                 inGal.quant75DensN=quantile(nDens,0.75);
            %                 inGal.stdDensN=std(nDens);
            
            
            inGal.velDisp(k,i)=calc_standardDev(vv);
            inGal.velDispMW(k,i)=calc_standardDev(vv,mm);
            
            
            inGal.gasMass(k,i)=sum(mm);
            inGal.avgDensN(k,i)=sum(mm)/sum(mm./nDens); % in cm^-3
            
            inGal.cellNum(k,i)=sum(mask);
            
            
            if fullSnapMask(snap+1)
                tc=tcool(mask);
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                
                
                inGal.meanTcMW(k,cntExt)=sum(mm.*tc)/sum(mm);
                inGal.medianTc(k,cntExt)=median(tc);
                %                 inGal.quant25Tc(k,i)=quantile(tc,0.25);
                %                 inGal.quant75Tc(k,i)=quantile(tc,0.75);
                %                 inGal.stdTc(k,i)=std(tc);
                inGal.EnergyDissipation(k,cntExt)=sum(enrgyDiss);
                inGal.meanMachEW(k,cntExt)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                inGal.meanMach(k,cntExt)=mean(mach(mach>0));
                inGal.medianMach(k,cntExt)=median(mach(mach>0));
                %                 inGal.quant25Mach(k,i)=quantile(mach(mach>0),0.25);
                %                 inGal.quant75Mach(k,i)=quantile(mach(mach>0),0.75);
                %                 inGal.stdMach(k,i)=std(mach(mach>0));
                
                
            end
            
            
            % mass in different gas components
            if k==2
                % star forming gas
                maskt=mask & sfMask;
                inGal.sfGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inGal.sfGas.EnergyDissipation(cntExt)=sum(ediss);
                    inGal.sfGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % cold & dense gas
                maskt=mask & coldDenseMask;
                inGal.coldDenseGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inGal.coldDenseGas.EnergyDissipation(cntExt)=sum(ediss);
                    inGal.coldDenseGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % cold & dilute gass
                maskt=mask & coldDiluteDMask;
                inGal.coldDiluteGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inGal.coldDiluteGas.EnergyDissipation(cntExt)=sum(ediss);
                    inGal.coldDiluteGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % warm & hot gass
                maskt=mask & warmHotMask;
                inGal.warmHotGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inGal.warmHotGas.EnergyDissipation(cntExt)=sum(ediss);
                    inGal.warmHotGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % hot gas
                maskt=mask & hotMask;
                inGal.hotGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inGal.hotGas.EnergyDissipation(cntExt)=sum(ediss);
                    inGal.hotGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                                               
            end
            
            
            
            
        end
        
    end
    
    
    
    %% for gas beyond on th galaxy and within rhalfGas
    distMask=gasDist>2.0.*rhalfStar &  ...
        gasDist<=rhalfGas;
    for k=1:2
        switch k
            case 1
                mask=distMask & ~sfMask;
            case 2
                mask=distMask;
        end
        
        
        
        if any(mask)
            mm=mass(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            vv=vel(mask);
            nDens=nDensity(mask);
            
            
            inCGM.meanTempMW(k,i)=sum(mm.*tmp)/sum(mm);
            inCGM.medianTemp(k,i)=median(tmp);
            
            inCGM.meanEntMW(k,i)=sum(mm.*ent)/sum(mm);
            inCGM.medianEnt(k,i)=median(ent);
            
            inCGM.meanDensN(k,i)=mean(nDens);
            inCGM.medianDensN(k,i)=median(nDens);
            
            inCGM.velDisp(k,i)=calc_standardDev(vv);
            inCGM.velDispMW(k,i)=calc_standardDev(vv,mm);
            
            
            inCGM.gasMass(k,i)=sum(mm);
            inCGM.avgDensN(k,i)=sum(mm)/sum(mm./nDens); % in cm^-3
            
            inCGM.cellNum(k,i)=sum(mask);
            
            
            if fullSnapMask(snap+1)
                tc=tcool(mask);
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                
                inCGM.meanTcMW(k,cntExt)=sum(mm.*tc)/sum(mm);
                inCGM.medianTc(k,cntExt)=median(tc);
                inCGM.EnergyDissipation(k,cntExt)=sum(enrgyDiss);
                inCGM.meanMachEW(k,cntExt)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                inCGM.meanMach(k,cntExt)=mean(mach(mach>0));
                inCGM.medianMach(k,cntExt)=median(mach(mach>0));
                
            end
            
            
            % mass in different gas components
            if k==2
                % star forming gas
                maskt=mask & sfMask;
                inCGM.sfGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inCGM.sfGas.EnergyDissipation(cntExt)=sum(ediss);
                    inCGM.sfGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % cold & dense gas
                maskt=mask & coldDenseMask;
                inCGM.coldDenseGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inCGM.coldDenseGas.EnergyDissipation(cntExt)=sum(ediss);
                    inCGM.coldDenseGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % cold & dilute gass
                maskt=mask & coldDiluteDMask;
                inCGM.coldDiluteGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inCGM.coldDiluteGas.EnergyDissipation(cntExt)=sum(ediss);
                    inCGM.coldDiluteGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % warm & hot gass
                maskt=mask & warmHotMask;
                inCGM.warmHotGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inCGM.warmHotGas.EnergyDissipation(cntExt)=sum(ediss);
                    inCGM.warmHotGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % hot gas
                maskt=mask & hotMask;
                inCGM.hotGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inCGM.hotGas.EnergyDissipation(cntExt)=sum(ediss);
                    inCGM.hotGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                
                
            end
            
            
            
            
        end
        
    end
    
    
    
    
    %% for all gas beyond on the galaxy
    distMask=gasDist>2.0.*rhalfStar ;
    
    
    for k=1:2
        switch k
            case 1
                mask=distMask & ~sfMask;
            case 2
                mask=distMask;
        end
        
        
        
        if any(mask)
            mm=mass(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            vv=vel(mask);
            nDens=nDensity(mask);
            
            inSub.meanTempMW(k,i)=sum(mm.*tmp)/sum(mm);
            inSub.medianTemp(k,i)=median(tmp);
            
            inSub.meanEntMW(k,i)=sum(mm.*ent)/sum(mm);
            inSub.medianEnt(k,i)=median(ent);
            
            inSub.meanDensN(k,i)=mean(nDens);
            inSub.medianDensN(k,i)=median(nDens);
            
            
            inSub.velDisp(k,i)=calc_standardDev(vv);
            inSub.velDispMW(k,i)=calc_standardDev(vv,mm);
            
            
            inSub.gasMass(k,i)=sum(mm);
            inSub.avgDensN(k,i)=sum(mm)/sum(mm./nDens); % in cm^-3
            
            inSub.cellNum(k,i)=sum(mask);
            
            
            if fullSnapMask(snap+1)
                tc=tcool(mask);
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                
                inSub.meanTcMW(k,cntExt)=sum(mm.*tc)/sum(mm);
                inSub.medianTc(k,cntExt)=median(tc);
                inSub.EnergyDissipation(k,cntExt)=sum(enrgyDiss);
                inSub.meanMachEW(k,cntExt)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                inSub.meanMach(k,cntExt)=mean(mach(mach>0));
                inSub.medianMach(k,cntExt)=median(mach(mach>0));
                
                
            end
            
            
            % mass in different gas components
            if k==2
                % star forming gas
                maskt=mask & sfMask;
                inSub.sfGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inSub.sfGas.EnergyDissipation(cntExt)=sum(ediss);
                    inSub.sfGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % cold & dense gas
                maskt=mask & coldDenseMask;
                inSub.coldDenseGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inSub.coldDenseGas.EnergyDissipation(cntExt)=sum(ediss);
                    inSub.coldDenseGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % cold & dilute gass
                maskt=mask & coldDiluteDMask;
                inSub.coldDiluteGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inSub.coldDiluteGas.EnergyDissipation(cntExt)=sum(ediss);
                    inSub.coldDiluteGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % warm & hot gass
                maskt=mask & warmHotMask;
                inSub.warmHotGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inSub.warmHotGas.EnergyDissipation(cntExt)=sum(ediss);
                    inSub.warmHotGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                
                % hot gas
                maskt=mask & hotMask;
                inSub.hotGas.mass(i)=sum(mass(maskt));
                
                if fullSnapMask(snap+1)
                    ediss=gas.EnergyDissipation(maskt);
                    inSub.hotGas.EnergyDissipation(cntExt)=sum(ediss);
                    inSub.hotGas.meanMachEW(cntExt)=sum(gas.Machnumber(maskt).*ediss)/sum(ediss);
                end
                          
                
            end
            
            
            
            
        end
        
    end
    
    
    %% BH stuff
    
    % load bh data
    bhs=illustris.snapshot.loadSubhalo(BASEPATH, tree.SnapNum(i), tree.SubfindID(i),'bh',...
        {'BH_CumEgyInjection_QM','BH_CumEgyInjection_RM','BH_Mass','BH_Mdot','Coordinates'});
    
    if bhs.count==0
        continue
    end
    
    
    bhs.newCoord = illustris.utils.centerObject(bhs.Coordinates,tree.SubhaloPos(:,i));
    bhsDist=sqrt( sum(bhs.newCoord.^2,1));
    
    
    %% for bhs within the gal (2*r_half - stellar)
    
    mask=bhsDist<=2.0.*rhalfStar;
    
    if any(mask)
        inGal.cumEngQM(i)=sum(bhs.BH_CumEgyInjection_QM(mask)).*illUnits.BHEnergyFactor;
        inGal.cumEngRM(i)=sum(bhs.BH_CumEgyInjection_RM(mask)).*illUnits.BHEnergyFactor;
        [mx, ix]=max(bhs.BH_Mass(mask));
        inGal.bhMassMax(i)=mx.*illUnits.massUnit;
        inGal.bhMassSum(i)=sum(bhs.BH_Mass(mask)).*illUnits.massUnit;
        inGal.bhMdot(i)=bhs.BH_Mdot(ix).*illUnits.BHMdotFactor;
        inGal.bhCount(i)=sum(mask);
        
    end
    
    
    %% for gas within the cgm (2*r_half  - r_half gas)
    
    mask=bhsDist>2.0.*rhalfStar &  ...
        bhsDist<=rhalfGas;
    
    if any(mask)
        inCGM.cumEngQM(i)=sum(bhs.BH_CumEgyInjection_QM(mask)).*illUnits.BHEnergyFactor;
        inCGM.cumEngRM(i)=sum(bhs.BH_CumEgyInjection_RM(mask)).*illUnits.BHEnergyFactor;
        [mx, ix]=max(bhs.BH_Mass(mask));
        inCGM.bhMassMax(i)=mx.*illUnits.massUnit;
        inCGM.bhMdot(i)=bhs.BH_Mdot(ix).*illUnits.BHMdotFactor;
        inCGM.bhMassSum(i)=sum(bhs.BH_Mass(mask)).*illUnits.massUnit;
        inCGM.bhCount(i)=sum(mask);
        
    end
    
    %% gas within 2* half stellar mass radius and edge
    mask=bhsDist>2.0.*rhalfStar ;
    
    if any(mask)
        inSub.cumEngQM(i)=sum(bhs.BH_CumEgyInjection_QM(mask)).*illUnits.BHEnergyFactor;
        inSub.cumEngRM(i)=sum(bhs.BH_CumEgyInjection_RM(mask)).*illUnits.BHEnergyFactor;
        [mx, ix]=max(bhs.BH_Mass(mask));
        inSub.bhMassMax(i)=mx.*illUnits.massUnit;
        inSub.bhMdot(i)=bhs.BH_Mdot(ix).*illUnits.BHMdotFactor;
        inSub.bhMassSum(i)=sum(bhs.BH_Mass(mask)).*illUnits.massUnit;
        inSub.bhCount(i)=sum(mask);
        
    end
    
end

res.inGal=inGal;
res.inCGM=inCGM;
res.inSub=inSub;




