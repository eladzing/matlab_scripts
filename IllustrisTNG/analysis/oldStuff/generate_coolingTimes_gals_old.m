%% Study the differences between the gas in sf vs. quiescent centals
% go over all galaxies and extract cooling times for the gas and entropy as
% well
%
% several areas: within 2*r_h, all gas in the subhalo, others?

%% set framework

%global matFilePath
%global cosmoStruct
global illUnits

sim='100';
snap=99; %z=0
bp=illustris.set_env(sim,'draco');


if ~exist('readFlag','var')
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

massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask= subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas & subs.SubhaloMassInRadType(1,:)>0;
tCoolStruct.galMass=massAllGals;
tCoolStruct.galMask=galMask;
%% generate birds
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;
for id=0:len-1
    
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if ~galMask(id+1)
        continue
    end
    
    cnt=cnt+1;
    % load gas from in sub halo
    gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas',{'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','GFM_CoolingRate'});
    
    
    
    % get density and temperature
    try 
        dens=gas.Density.*illUnits.numberDensityFactor; %in cm^-3
    catch
        fprintf('subhalo id: %s \n',num2str(id))
    end
       
    gas=illustris.utils.addTemperature( gas );
    %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
    
    
    mass=gas.Masses.*illUnits.massUnit; %in Solarmass
    tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr^-1
    gas=illustris.utils.addEntropy( gas );
    
    
    % find distance from galaxy center
    
    gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
    gasDist=sqrt( sum(gas.newCoord.^2,1));
    
    
    % find importat radii
    rmax=max(gasDist);
    rhalfStar=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1); % stellar half mass radius
    rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1); % gas half mass radius
    tCoolStruct.rmax(cnt)=rmax.*illUnits.lengthUnit;
    tCoolStruct.rhalfStar(cnt)=rhalfStar.*illUnits.lengthUnit;% in kpc
    tCoolStruct.rhalfGas(cnt)=rhalfGas.*illUnits.lengthUnit; % in kpc
    
    % identify Star-forming gas
    phaseDiagram_polygons
    sfMask=inpolygon(log10(dens),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
    
    
    %% for gas within the gal (2*r_half - stellar)
    
    distMask=gasDist<=2.0.*rhalfStar;
    
    for i=1:2
        switch i
            case 1
                mask=distMask & ~sfMask;
            case 2
                mask=distMask;
        end
        
        
        
        if any(mask)
            mm=mass(mask);
            tc=tcool(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            
            
            tCoolStruct.inGal.meanTc(cnt,i)=mean(tc);
            tCoolStruct.inGal.meanTcMW(cnt,i)=sum(mm.*tc)/sum(mm);
            tCoolStruct.inGal.quantileTc(1:5,cnt,i)=quantile(tc,[0.1 0.25 0.5 0.75 0.9]);
            tCoolStruct.inGal.stdTc(cnt,i)=std(tc);
            
            tCoolStruct.inGal.meanTemp(cnt,i)=mean(tmp);
            tCoolStruct.inGal.meanTempMW(cnt,i)=sum(mm.*tmp)/sum(mm);
            tCoolStruct.inGal.mass(cnt,i)=sum(mm);
            tCoolStruct.inGal.num(cnt,i)=sum(mask);
            
            tCoolStruct.inGal.meanEnt(cnt,i)=mean(ent);
            tCoolStruct.inGal.meanEntMW(cnt,i)=sum(mm.*ent)/sum(mm);
            tCoolStruct.inGal.quantileEnt(1:5,cnt,i)=quantile(ent,[0.1 0.25 0.5 0.75 0.9]);
        else
            
            tCoolStruct.inGal.meanTc(cnt,i)=0;
            tCoolStruct.inGal.meanTcMW(cnt,i)=0;
            tCoolStruct.inGal.quantileTc(1:5,cnt,i)=0;
            tCoolStruct.inGal.stdTc(cnt,i)=0;
            
            tCoolStruct.inGal.meanTemp(cnt,i)=0;
            tCoolStruct.inGal.meanTempMW(cnt,i)=0;
            tCoolStruct.inGal.mass(cnt,i)=0;
            tCoolStruct.inGal.num(cnt,i)=0;
            
            tCoolStruct.inGal.meanEnt(cnt,i)=0;
            tCoolStruct.inGal.meanEntMW(cnt,i)=0;
            tCoolStruct.inGal.quantileEnt(1:5,cnt,i)=0;
        end
        
    end
    
    
    
    %% for gas within the cgm (2*r_half  - r_half gas)
    
    distMask=gasDist>2.0.*rhalfStar &  ...
        gasDist<=rhalfGas;
    
    for i=1:2
        switch i
            case 1
                mask=distMask & ~sfMask;
            case 2
                mask=distMask;
        end
        
        
        if any(mask)
            mm=mass(mask);
            tc=tcool(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            
            tCoolStruct.inCGM.meanTc(cnt,i)=mean(tc);
            tCoolStruct.inCGM.meanTcMW(cnt,i)=sum(mm.*tc)/sum(mm);
            tCoolStruct.inCGM.quantileTc(1:5,cnt,i)=quantile(tc,[0.1 0.25 0.5 0.75 0.9]);
            tCoolStruct.inCGM.stdTc(cnt,i)=std(tc);
            
            
            tCoolStruct.inCGM.meanTemp(cnt,i)=mean(tmp);
            tCoolStruct.inCGM.meanTempMW(cnt,i)=sum(mm.*tmp)/sum(mm);
            tCoolStruct.inCGM.mass(cnt,i)=sum(mm);
            tCoolStruct.inCGM.num(cnt,i)=sum(mask);
            
            tCoolStruct.inCGM.meanEnt(cnt,i)=mean(ent);
            tCoolStruct.inCGM.meanEntMW(cnt,i)=sum(mm.*ent)/sum(mm);
            tCoolStruct.inCGM.quantileEnt(1:5,cnt,i)=quantile(ent,[0.1 0.25 0.5 0.75 0.9]);
            
        else
            tCoolStruct.inCGM.meanTc(cnt,i)=0;
            tCoolStruct.inCGM.meanTcMW(cnt,i)=0;
            tCoolStruct.inCGM.quantileTc(1:5,cnt,i)=0;
            tCoolStruct.inCGM.stdTc(cnt,i)=0;
            
            
            tCoolStruct.inCGM.meanTemp(cnt,i)=0;
            tCoolStruct.inCGM.meanTempMW(cnt,i)=0;
            tCoolStruct.inCGM.mass(cnt,i)=0;
            tCoolStruct.inCGM.num(cnt,i)=0;
            
            tCoolStruct.inCGM.meanEnt(cnt,i)=0;
            tCoolStruct.inCGM.meanEntMW(cnt,i)=0;
            tCoolStruct.inCGM.quantileEnt(1:5,cnt,i)=0;
            
        end
    end
    
    
    %% gas within 2* half stellar mass radius and edge
    distMask=gasDist>2.0.*rhalfStar ;
    for i=1:2
        switch i
            case 1
                mask=distMask & ~sfMask;
            case 2
                mask=distMask;
        end
        
        if any(mask)
            
            mm=mass(mask);
            tc=tcool(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            
            tCoolStruct.inSub.meanTc(cnt,i)=mean(tc);
            tCoolStruct.inSub.meanTcMW(cnt,i)=sum(mm.*tc)/sum(mm);
            tCoolStruct.inSub.quantileTc(1:5,cnt,i)=quantile(tc,[0.1 0.25 0.5 0.75 0.9]);
            tCoolStruct.inSub.stdTc(cnt,i)=std(tc);
            
            tCoolStruct.inSub.meanTemp(cnt,i)=mean(tmp);
            tCoolStruct.inSub.meanTempMW(cnt,i)=sum(mm.*tmp)/sum(mm);
            tCoolStruct.inSub.mass(cnt,i)=sum(mm);
            tCoolStruct.inSub.num(cnt,i)=sum(mask);
            
            tCoolStruct.inSub.meanEnt(cnt,i)=mean(ent);
            tCoolStruct.inSub.meanEntMW(cnt,i)=sum(mm.*ent)/sum(mm);
            tCoolStruct.inSub.quantileEnt(1:5,cnt,i)=quantile(ent,[0.1 0.25 0.5 0.75 0.9]);
            
        else
            tCoolStruct.inSub.meanTc(cnt,i)=0;
            tCoolStruct.inSub.meanTcMW(cnt,i)=0;
            tCoolStruct.inSub.quantileTc(1:5,cnt,i)=0;
            tCoolStruct.inSub.stdTc(cnt,i)=0;
            
            tCoolStruct.inSub.meanTemp(cnt,i)=0;
            tCoolStruct.inSub.meanTempMW(cnt,i)=0;
            tCoolStruct.inSub.mass(cnt,i)=0;
            tCoolStruct.inSub.num(cnt,i)=0;
            
            tCoolStruct.inSub.meanEnt(cnt,i)=0;
            tCoolStruct.inSub.meanEntMW(cnt,i)=0;
            tCoolStruct.inSub.quantileEnt(1:5,cnt,i)=0;
            
        end
        
        
    end
    
    
    
    
    
    
end


%% save to mat file
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    fname=sprintf('cooling_times_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'tCoolStruct','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
end
fprintf(' *** DONE!  *** \n');


