%% Study the differences between the gas in sf vs. quiescent centals
% go over all galaxies and extract cooling times for the gas and entropy as
% well
%
% several areas: within 2*r_h, all gas in the subhalo, others?

%% set framework

%global matFilePath
%global cosmoStruct


%sim='100';
snap=99; %z=0
bp=illustris.set_env(simName); %,'draco');

global illUnits

if ~exist('readFlag','var')
    readFlag=true;
end

if readFlag
    fprintf(' *** Reading data *** \n');
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
end


massThresh=10^7; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
units; % load general unit structure in cgs.
%% select galaxies

%massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
%massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
%galMask= subs.SubhaloFlag & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh);
%galMask=galMask & subsInfo.hasGas;
%tCoolStruct.galMass=massAllGals;
%tCoolStruct.galMask=galMask;

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;




%% initialize to zero

tCoolStruct.inGal.meanTc=zeros(2,len);
tCoolStruct.inGal.meanTcMW=zeros(2,len);
tCoolStruct.inGal.medianTc=zeros(2,len);
tCoolStruct.inGal.quant25Tc=zeros(2,len);
tCoolStruct.inGal.quant75Tc=zeros(2,len);
tCoolStruct.inGal.stdTc=zeros(2,len);

tCoolStruct.inGal.meanTemp=zeros(2,len);
tCoolStruct.inGal.meanTempMW=zeros(2,len);
tCoolStruct.inGal.medianTemp=zeros(2,len);
tCoolStruct.inGal.quant25Temp=zeros(2,len);
tCoolStruct.inGal.quant75Temp=zeros(2,len);
tCoolStruct.inGal.stdTemp=zeros(2,len);

tCoolStruct.inGal.meanEnt=zeros(2,len);
tCoolStruct.inGal.meanEntMW=zeros(2,len);
tCoolStruct.inGal.medianEnt=zeros(2,len);
tCoolStruct.inGal.quant25Ent=zeros(2,len);
tCoolStruct.inGal.quant75Ent=zeros(2,len);
tCoolStruct.inGal.stdEnt=zeros(2,len);

tCoolStruct.inGal.meanDensN=zeros(2,len);
tCoolStruct.inGal.medianDensN=zeros(2,len);
tCoolStruct.inGal.quant25DensN=zeros(2,len);
tCoolStruct.inGal.quant75DensN=zeros(2,len);
tCoolStruct.inGal.stdDensN=zeros(2,len);

tCoolStruct.inGal.EnergyDissipation=zeros(2,len);
tCoolStruct.inGal.meanMachEW=zeros(2,len);
tCoolStruct.inGal.meanMach=zeros(2,len);
tCoolStruct.inGal.medianMach=zeros(2,len);
tCoolStruct.inGal.quant25Mach=zeros(2,len);
tCoolStruct.inGal.quant75Mach=zeros(2,len);
tCoolStruct.inGal.stdMach=zeros(2,len);
tCoolStruct.inGal.velDisp=zeros(2,len);
tCoolStruct.inGal.velDispMW=zeros(2,len);

tCoolStruct.inGal.gasMass=zeros(2,len);
tCoolStruct.inGal.avgDensN=zeros(2,len);
tCoolStruct.inGal.cellNum=zeros(2,len);

tCoolStruct.inCGM=tCoolStruct.inGal;
tCoolStruct.inSub=tCoolStruct.inGal;
tCoolStruct.inOut=tCoolStruct.inGal;

for id=0:len-1
    
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if galMask(id+1)
        
        
        %cnt=cnt+1;
        % load gas from in sub halo
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas',{'Coordinates','Density','ElectronAbundance',...
            'InternalEnergy','Masses','GFM_CoolingRate','EnergyDissipation','Machnumber','Velocities',...
            'StarFormationRate'});
        
        if gas.count==0
            continue
        end
        
        % get density and temperature
        try
            nDensity=gas.Density.*illUnits.numberDensityFactor; %in cm^-3
        catch
            fprintf('subhalo id: %s \n',num2str(id))
        end
        
        gas=illustris.utils.addTemperature( gas );
        %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
        
        
        mass=gas.Masses.*illUnits.massUnit; %in Solarmass
        tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr^-1
        gas=illustris.utils.addEntropy( gas );
        
        % deal with velocitites
        
        for k=1:3
            vv=gas.Velocities(k,:).*sqrt(illustris.utils.snap2aexpn(snap))-subs.SubhaloVel(k,id+1);
        end
        vel=sqrt(sum(vv.^2,1));
        clear vv
        
        
        
        % find distance from galaxy center
        
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
        gasDist=sqrt( sum(gas.newCoord.^2,1));
        
        
        % find importat radii
        rmax=max(gasDist);
        rhalfStar=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1); % stellar half mass radius
        rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1); % gas half mass radius
        
        %tCoolStruct.rmax(id+)=rmax.*illUnits.lengthUnit;
        %tCoolStruct.rhalfStar(cnt)=rhalfStar.*illUnits.lengthUnit;% in kpc
        %tCoolStruct.rhalfGas(cnt)=rhalfGas.*illUnits.lengthUnit; % in kpc
        
        % identify Star-forming gas
        polys=phaseDiagram_polygons('noplot');
        sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
        
        
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
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                vv=vel(mask);
                nDens=nDensity(mask);
                
                tCoolStruct.inGal.meanTc(i,id+1)=mean(tc);
                tCoolStruct.inGal.meanTcMW(i,id+1)=sum(mm.*tc)/sum(mm);
                tCoolStruct.inGal.medianTc(i,id+1)=median(tc);
                tCoolStruct.inGal.quant25Tc(i,id+1)=quantile(tc,0.25);
                tCoolStruct.inGal.quant75Tc(i,id+1)=quantile(tc,0.75);
                tCoolStruct.inGal.stdTc(i,id+1)=std(tc);
                
                tCoolStruct.inGal.meanTemp(i,id+1)=mean(tmp);
                tCoolStruct.inGal.meanTempMW(i,id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inGal.medianTemp(i,id+1)=median(tmp);
                tCoolStruct.inGal.quant25Temp(i,id+1)=quantile(tmp,0.25);
                tCoolStruct.inGal.quant75Temp(i,id+1)=quantile(tmp,0.75);
                tCoolStruct.inGal.stdTemp(i,id+1)=std(tmp);
                
                tCoolStruct.inGal.meanEnt(i,id+1)=mean(ent);
                tCoolStruct.inGal.meanEntMW(i,id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inGal.medianEnt(i,id+1)=median(ent);
                tCoolStruct.inGal.quant25Ent(i,id+1)=quantile(ent,0.25);
                tCoolStruct.inGal.quant75Ent(i,id+1)=quantile(ent,0.75);
                tCoolStruct.inGal.stdEnt(i,id+1)=std(ent);
                
                tCoolStruct.inGal.meanDensN(i,id+1) =mean(nDens);
                tCoolStruct.inGal.medianDensN(i,id+1) =median(nDens);
                tCoolStruct.inGal.quant25DensN(i,id+1) =quantile(nDens,0.25);
                tCoolStruct.inGal.quant75DensN(i,id+1) =quantile(nDens,0.75);
                tCoolStruct.inGal.stdDensN(i,id+1) =std(nDens);
                
                tCoolStruct.inGal.EnergyDissipation(i,id+1)=sum(enrgyDiss);
                tCoolStruct.inGal.meanMachEW(i,id+1)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                tCoolStruct.inGal.meanMach(i,id+1)=mean(mach(mach>0));
                tCoolStruct.inGal.medianMach(i,id+1)=median(mach(mach>0));
                tCoolStruct.inGal.quant25Mach(i,id+1)=quantile(mach(mach>0),0.25);
                tCoolStruct.inGal.quant75Mach(i,id+1)=quantile(mach(mach>0),0.75);
                tCoolStruct.inGal.stdMach(i,id+1)=std(mach(mach>0));
                tCoolStruct.inGal.velDisp(i,id+1)=calc_standardDev(vv);
                tCoolStruct.inGal.velDispMW(i,id+1)=calc_standardDev(vv,mm);
                
                
                tCoolStruct.inGal.gasMass(i,id+1)=sum(mm);
                tCoolStruct.inGal.avgDensN(i,id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                
                tCoolStruct.inGal.cellNum(i,id+1)=sum(mask);
                
                
                %            tCoolStruct.inGal.quantileEnt(1:5,cnt,i)=quantile(ent,[0.1 0.25 0.5 0.75 0.9]);
                
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
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                vv=vel(mask);
                nDens=nDensity(mask);
                
                tCoolStruct.inCGM.meanTc(i,id+1)=mean(tc);
                tCoolStruct.inCGM.meanTcMW(i,id+1)=sum(mm.*tc)/sum(mm);
                tCoolStruct.inCGM.medianTc(i,id+1)=median(tc);
                tCoolStruct.inCGM.quant25Tc(i,id+1)=quantile(tc,0.25);
                tCoolStruct.inCGM.quant75Tc(i,id+1)=quantile(tc,0.75);
                tCoolStruct.inCGM.stdTc(i,id+1)=std(tc);
                
                tCoolStruct.inCGM.meanTemp(i,id+1)=mean(tmp);
                tCoolStruct.inCGM.meanTempMW(i,id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inCGM.medianTemp(i,id+1)=median(tmp);
                tCoolStruct.inCGM.quant25Temp(i,id+1)=quantile(tmp,0.25);
                tCoolStruct.inCGM.quant75Temp(i,id+1)=quantile(tmp,0.75);
                tCoolStruct.inCGM.stdTemp(i,id+1)=std(tmp);
                
                tCoolStruct.inCGM.meanEnt(i,id+1)=mean(ent);
                tCoolStruct.inCGM.meanEntMW(i,id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inCGM.medianEnt(i,id+1)=median(ent);
                tCoolStruct.inCGM.quant25Ent(i,id+1)=quantile(ent,0.25);
                tCoolStruct.inCGM.quant75Ent(i,id+1)=quantile(ent,0.75);
                tCoolStruct.inCGM.stdEnt(i,id+1)=std(ent);
                
                tCoolStruct.inCGM.meanDensN(i,id+1) =mean(nDens);
                tCoolStruct.inCGM.medianDensN(i,id+1) =median(nDens);
                tCoolStruct.inCGM.quant25DensN(i,id+1) =quantile(nDens,0.25);
                tCoolStruct.inCGM.quant75DensN(i,id+1) =quantile(nDens,0.75);
                tCoolStruct.inCGM.stdDensN(i,id+1) =std(nDens);
                
                
                tCoolStruct.inCGM.EnergyDissipation(i,id+1)=sum(enrgyDiss);
                tCoolStruct.inCGM.meanMachEW(i,id+1)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                tCoolStruct.inCGM.meanMach(i,id+1)=mean(mach(mach>0));
                tCoolStruct.inCGM.medianMach(i,id+1)=median(mach(mach>0));
                tCoolStruct.inCGM.quant25Mach(i,id+1)=quantile(mach(mach>0),0.25);
                tCoolStruct.inCGM.quant75Mach(i,id+1)=quantile(mach(mach>0),0.75);
                tCoolStruct.inCGM.stdMach(i,id+1)=std(mach(mach>0));
                tCoolStruct.inCGM.velDisp(i,id+1)=calc_standardDev(vv);
                tCoolStruct.inCGM.velDispMW(i,id+1)=calc_standardDev(vv,mm);
                
                tCoolStruct.inCGM.gasMass(i,id+1)=sum(mm);
                tCoolStruct.inCGM.avgDensN(i,id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                tCoolStruct.inCGM.cellNum(i,id+1)=sum(mask);
                
                
                
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
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                vv=vel(mask);
                nDens=nDensity(mask);
                
                
                tCoolStruct.inSub.meanTc(i,id+1)=mean(tc);
                tCoolStruct.inSub.meanTcMW(i,id+1)=sum(mm.*tc)/sum(mm);
                tCoolStruct.inSub.medianTc(i,id+1)=median(tc);
                tCoolStruct.inSub.quant25Tc(i,id+1)=quantile(tc,0.25);
                tCoolStruct.inSub.quant75Tc(i,id+1)=quantile(tc,0.75);
                tCoolStruct.inSub.stdTc(i,id+1)=std(tc);
                
                tCoolStruct.inSub.meanTemp(i,id+1)=mean(tmp);
                tCoolStruct.inSub.meanTempMW(i,id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inSub.medianTemp(i,id+1)=median(tmp);
                tCoolStruct.inSub.quant25Temp(i,id+1)=quantile(tmp,0.25);
                tCoolStruct.inSub.quant75Temp(i,id+1)=quantile(tmp,0.75);
                tCoolStruct.inSub.stdTemp(i,id+1)=std(tmp);
                
                tCoolStruct.inSub.meanEnt(i,id+1)=mean(ent);
                tCoolStruct.inSub.meanEntMW(i,id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inSub.medianEnt(i,id+1)=median(ent);
                tCoolStruct.inSub.quant25Ent(i,id+1)=quantile(ent,0.25);
                tCoolStruct.inSub.quant75Ent(i,id+1)=quantile(ent,0.75);
                tCoolStruct.inSub.stdEnt(i,id+1)=std(ent);
                
                tCoolStruct.inSub.meanDensN(i,id+1) =mean(nDens);
                tCoolStruct.inSub.medianDensN(i,id+1) =median(nDens);
                tCoolStruct.inSub.quant25DensN(i,id+1) =quantile(nDens,0.25);
                tCoolStruct.inSub.quant75DensN(i,id+1) =quantile(nDens,0.75);
                tCoolStruct.inSub.stdDensN(i,id+1) =std(nDens);
                
                tCoolStruct.inSub.EnergyDissipation(i,id+1)=sum(enrgyDiss);
                tCoolStruct.inSub.meanMachEW(i,id+1)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                tCoolStruct.inSub.meanMach(i,id+1)=mean(mach(mach>0));
                tCoolStruct.inSub.medianMach(i,id+1)=median(mach(mach>0));
                tCoolStruct.inSub.quant25Mach(i,id+1)=quantile(mach(mach>0),0.25);
                tCoolStruct.inSub.quant75Mach(i,id+1)=quantile(mach(mach>0),0.75);
                tCoolStruct.inSub.stdMach(i,id+1)=std(mach(mach>0));
                tCoolStruct.inSub.velDisp(i,id+1)=calc_standardDev(vv);
                tCoolStruct.inSub.velDispMW(i,id+1)=calc_standardDev(vv,mm);
                
                
                tCoolStruct.inSub.gasMass(i,id+1)=sum(mm);
                tCoolStruct.inSub.avgDensN(i,id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                tCoolStruct.inSub.cellNum(i,id+1)=sum(mask);
                
                
            end
            
            
        end
        
        
        %% gas from  half gas mass radius and edge
        distMask=gasDist>rhalfGas ;
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
                enrgyDiss=gas.EnergyDissipation(mask);
                mach=gas.Machnumber(mask);
                vv=vel(mask);
                nDens=nDensity(mask);
                
                
                tCoolStruct.inOut.meanTc(i,id+1)=mean(tc);
                tCoolStruct.inOut.meanTcMW(i,id+1)=sum(mm.*tc)/sum(mm);
                tCoolStruct.inOut.medianTc(i,id+1)=median(tc);
                tCoolStruct.inOut.quant25Tc(i,id+1)=quantile(tc,0.25);
                tCoolStruct.inOut.quant75Tc(i,id+1)=quantile(tc,0.75);
                tCoolStruct.inOut.stdTc(i,id+1)=std(tc);
                
                tCoolStruct.inOut.meanTemp(i,id+1)=mean(tmp);
                tCoolStruct.inOut.meanTempMW(i,id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inOut.medianTemp(i,id+1)=median(tmp);
                tCoolStruct.inOut.quant25Temp(i,id+1)=quantile(tmp,0.25);
                tCoolStruct.inOut.quant75Temp(i,id+1)=quantile(tmp,0.75);
                tCoolStruct.inOut.stdTemp(i,id+1)=std(tmp);
                
                tCoolStruct.inOut.meanEnt(i,id+1)=mean(ent);
                tCoolStruct.inOut.meanEntMW(i,id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inOut.medianEnt(i,id+1)=median(ent);
                tCoolStruct.inOut.quant25Ent(i,id+1)=quantile(ent,0.25);
                tCoolStruct.inOut.quant75Ent(i,id+1)=quantile(ent,0.75);
                tCoolStruct.inOut.stdEnt(i,id+1)=std(ent);
                
                tCoolStruct.inOut.meanDensN(i,id+1) =mean(nDens);
                tCoolStruct.inOut.medianDensN(i,id+1) =median(nDens);
                tCoolStruct.inOut.quant25DensN(i,id+1) =quantile(nDens,0.25);
                tCoolStruct.inOut.quant75DensN(i,id+1) =quantile(nDens,0.75);
                tCoolStruct.inOut.stdDensN(i,id+1) =std(nDens);
                
                tCoolStruct.inOut.EnergyDissipation(i,id+1)=sum(enrgyDiss);
                tCoolStruct.inOut.meanMachEW(i,id+1)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                tCoolStruct.inOut.meanMach(i,id+1)=mean(mach(mach>0));
                tCoolStruct.inOut.medianMach(i,id+1)=median(mach(mach>0));
                tCoolStruct.inOut.quant25Mach(i,id+1)=quantile(mach(mach>0),0.25);
                tCoolStruct.inOut.quant75Mach(i,id+1)=quantile(mach(mach>0),0.75);
                tCoolStruct.inOut.stdMach(i,id+1)=std(mach(mach>0));
                tCoolStruct.inOut.velDisp(i,id+1)=calc_standardDev(vv);
                tCoolStruct.inOut.velDispMW(i,id+1)=calc_standardDev(vv,mm);
                
                
                tCoolStruct.inOut.gasMass(i,id+1)=sum(mm);
                tCoolStruct.inOut.avgDensN(i,id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                tCoolStruct.inOut.cellNum(i,id+1)=sum(mask);
                
                
            end
            
            
        end
        
        
    end
    
    
    
    
    
end
tCoolStruct.galMask=galMask;

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

%% save to postprocessed catalog

%global myPOSTPROCESSING
% tCoolStruct.inGal.mask=floor(galMask);
% tCoolStruct.inCGM.mask=floor(galMask);
% tCoolStruct.inSub.mask=floor(galMask);
% illustris.utils.write_catalog( tCoolStruct.inGal,99,'name','gasProps_inGal','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( tCoolStruct.inCGM,99,'name','gasProps_inCGM','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( tCoolStruct.inSub,99,'name','gasProps_inSub','folderName','gasProperties','verbose');
%



