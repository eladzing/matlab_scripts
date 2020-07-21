%% this is an updated version of generate_coolingTimes_gals
% new things - medians and quantiles are now mass-median.
%            - STD is now mass weighted
%            - added mass distribution to some parameters
%            - commented out some values which have not been useful so far

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

illustris.utils.set_illUnits(snap);

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


massThresh=10^9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
units; % load general unit structure in cgs.
%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
%massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
%galMask= subs.SubhaloFlag & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh);
%galMask=galMask & subsInfo.hasGas;
tCoolStruct.galMass=massAllGals;
tCoolStruct.galMask=galMask;
massHist.galMass=massAllGals;
massHist.galMask=galMask;

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;

distLen=200;


%% initialize to zero

% tCoolStruct.inGal.meanTc=zeros(1,len);
tCoolStruct.inGal.meanTcMW=zeros(1,len);
tCoolStruct.inGal.medianTc=zeros(1,len);
tCoolStruct.inGal.quantTc=zeros(4,len);
%tCoolStruct.inGal.quant75Tc=zeros(1,len);
% tCoolStruct.inGal.stdTc=zeros(1,len);
tCoolStruct.inGal.stdTcMW=zeros(1,len);
massHist.inGal.tcool=zeros(distLen,len);
massHist.inGal.tcoolMass=zeros(distLen,len);

% tCoolStruct.inGal.meanTemp=zeros(1,len);
tCoolStruct.inGal.meanTempMW=zeros(1,len);
tCoolStruct.inGal.medianTemp=zeros(1,len);
tCoolStruct.inGal.quantTemp=zeros(4,len);
%tCoolStruct.inGal.quant75Temp=zeros(1,len);
% tCoolStruct.inGal.stdTemp=zeros(1,len);
tCoolStruct.inGal.stdTempMW=zeros(1,len);
massHist.inGal.temp=zeros(distLen,len);
massHist.inGal.tempMass=zeros(distLen,len);

% tCoolStruct.inGal.meanEnt=zeros(1,len);
tCoolStruct.inGal.meanEntMW=zeros(1,len);
tCoolStruct.inGal.medianEnt=zeros(1,len);
tCoolStruct.inGal.quantEnt=zeros(4,len);
%tCoolStruct.inGal.quant75Ent=zeros(1,len);
% tCoolStruct.inGal.stdEnt=zeros(1,len);
tCoolStruct.inGal.stdEntMW=zeros(1,len);
massHist.inGal.ent=zeros(distLen,len);
massHist.inGal.entMass=zeros(distLen,len);

tCoolStruct.inGal.meanDensN=zeros(1,len);
tCoolStruct.inGal.medianDensN=zeros(1,len);
tCoolStruct.inGal.quantDensN=zeros(4,len);
%tCoolStruct.inGal.quant75DensN=zeros(1,len);
tCoolStruct.inGal.stdDensN=zeros(1,len);
massHist.inGal.dens=zeros(distLen,len);
massHist.inGal.densMass=zeros(distLen,len);

% tCoolStruct.inGal.EnergyDissipation=zeros(2,len);
% tCoolStruct.inGal.meanMachEW=zeros(2,len);
% tCoolStruct.inGal.meanMach=zeros(2,len);
% tCoolStruct.inGal.medianMach=zeros(2,len);
% tCoolStruct.inGal.quant25Mach=zeros(2,len);
% tCoolStruct.inGal.quant75Mach=zeros(2,len);
% tCoolStruct.inGal.stdMach=zeros(2,len);
% tCoolStruct.inGal.velDisp=zeros(2,len);
% tCoolStruct.inGal.velDispMW=zeros(2,len);

tCoolStruct.inGal.gasMass=zeros(1,len);
tCoolStruct.inGal.sfrMass=zeros(1,len);
tCoolStruct.inGal.avgDensN=zeros(1,len);
tCoolStruct.inGal.cellNum=zeros(1,len);
massHist.inGal.gasMass=zeros(1,len);


tCoolStruct.inCGM=tCoolStruct.inGal;
tCoolStruct.inSub=tCoolStruct.inGal;
tCoolStruct.inOut=tCoolStruct.inGal;
massHist.inCGM=massHist.inGal;
massHist.inSub=massHist.inGal;
massHist.inOut=massHist.inGal;

qus=[0.1 0.25 0.5 0.75 0.9];

    
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
                'InternalEnergy','Masses','GFM_CoolingRate','StarFormationRate'});%'EnergyDissipation','Machnumber','Velocities',...
            
            
            if gas.count==0
                continue
            end
            
            % get density and temperature
            try
                nDensity=double(gas.Density.*illUnits.numberDensityFactor); %in cm^-3
            catch
                
                fprintf('subhalo id: %s \n',num2str(id))
            end
            
            gas=illustris.utils.addTemperature( gas );
            %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
            
            
            mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
            tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr^-1
            gas=illustris.utils.addEntropy( gas );
            
            % deal with velocitites
            
            %         for k=1:3
            %             vv=gas.Velocities(k,:).*sqrt(illustris.utils.snap2aexpn(snap))-subs.SubhaloVel(k,id+1);
            %         end
            %         vel=sqrt(sum(vv.^2,1));
            %         clear vv
            
            
            
            % find distance from galaxy center
            
            gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
            gasDist=sqrt( sum(double(gas.newCoord).^2,1));
            
            
            % find importat radii
            rmax=max(gasDist);
            rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)); % stellar half mass radius
            rhalfGas=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1)); % gas half mass radius
            
            %tCoolStruct.rmax(id+)=rmax.*illUnits.lengthUnit;
            %tCoolStruct.rhalfStar(cnt)=rhalfStar.*illUnits.lengthUnit;% in kpc
            %tCoolStruct.rhalfGas(cnt)=rhalfGas.*illUnits.lengthUnit; % in kpc
            
            % identify Star-forming gas
            %polys=phaseDiagram_polygons('noplot');
            sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
            
            
            %% for gas within the gal (2*r_half - stellar)
            
            distMask=gasDist<=2.0.*rhalfStar;
            mask=distMask & ~sfMask;
            
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);
                tmp=gas.Temperature(mask);
                ent=gas.Entropy(mask);
                %             enrgyDiss=gas.EnergyDissipation(mask);
                %             mach=gas.Machnumber(mask);
                %             vv=vel(mask);
                nDens=nDensity(mask);
                
                
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                if sum(tcMask)>0
                    [xx, massDist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus);
                    mtemp=massDist.*sum(mm(tcMask));
                    mtemp(end)=mtemp(end)+sum(mm(~tcMask));
                    massDist=mtemp./sum(mtemp);
                    
                    massHist.inGal.tcool(:,id+1)=xx;
                    massHist.inGal.tcoolMass(:,id+1)=massDist;
                    
                    tCoolStruct.inGal.meanTcMW(id+1)=sum(mm.*tc)/sum(mm);
                    tCoolStruct.inGal.stdTcMW(id+1)=calc_standardDev(tmp,mm);
                    tCoolStruct.inGal.medianTc(id+1)=10.^mus(3);
                    tCoolStruct.inGal.quantTc(:,id+1)=10.^mus([1 2 4 5]);
                end
                
                % temp
                [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
                
                massHist.inGal.temp(:,id+1)=xx;
                massHist.inGal.tempMass(:,id+1)=massDist;
                
                tCoolStruct.inGal.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inGal.stdTempMW(id+1)=calc_standardDev(tmp,mm);
                tCoolStruct.inGal.medianTemp(id+1)=10.^mus(3);
                tCoolStruct.inGal.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
                
                %entropy
                [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
                
                massHist.inGal.ent(:,id+1)=xx;
                massHist.inGal.entMass(:,id+1)=massDist;
                
                tCoolStruct.inGal.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inGal.stdEntMW(id+1)=calc_standardDev(ent,mm);
                tCoolStruct.inGal.medianEnt(id+1)=10.^mus(3);
                tCoolStruct.inGal.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
                
                %number density
                [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
                
                massHist.inGal.dens(:,id+1)=xx;
                massHist.inGal.densMass(:,id+1)=massDist;
                
                tCoolStruct.inGal.meanDensN(id+1)=mean(nDens);
                tCoolStruct.inGal.stdDensN(id+1)=std(nDens);
                tCoolStruct.inGal.medianDensN(id+1)=10.^mus(3);
                tCoolStruct.inGal.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
                
                
                
                %                 tCoolStruct.inGal.EnergyDissipation(i,id+1)=sum(enrgyDiss);
                %                 tCoolStruct.inGal.meanMachEW(i,id+1)=sum(mach.*enrgyDiss)/sum(enrgyDiss);
                %                 tCoolStruct.inGal.meanMach(i,id+1)=mean(mach(mach>0));
                %                 tCoolStruct.inGal.medianMach(i,id+1)=median(mach(mach>0));
                %                 tCoolStruct.inGal.quant25Mach(i,id+1)=quantile(mach(mach>0),0.25);
                %                 tCoolStruct.inGal.quant75Mach(i,id+1)=quantile(mach(mach>0),0.75);
                %                 tCoolStruct.inGal.stdMach(i,id+1)=std(mach(mach>0));
                %                 tCoolStruct.inGal.velDisp(i,id+1)=calc_standardDev(vv);
                %                 tCoolStruct.inGal.velDispMW(i,id+1)=calc_standardDev(vv,mm);
                
                
                tCoolStruct.inGal.gasMass(id+1)=sum(mm);
                tCoolStruct.inGal.sfrMass(id+1)=sum(mass(distMask & sfMask));
                tCoolStruct.inGal.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                
                tCoolStruct.inGal.cellNum(id+1)=sum(mask);
                
                massHist.inGal.gasMass(id+1)=sum(mm);
                
                
            end
            
            
            
            %% for gas within the cgm (2*r_half  - r_half gas)
            
            distMask=gasDist>2.0.*rhalfStar &  ...
                gasDist<=rhalfGas;
            
            mask=distMask & ~sfMask;
            
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);
                tmp=gas.Temperature(mask);
                ent=gas.Entropy(mask);
                nDens=nDensity(mask);
                
                
                %tcool
                
                tcMask=tc>0;
                if sum(tcMask)>0
                    [xx, massDist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus);
                    mtemp=massDist.*sum(mm(tcMask));
                    mtemp(end)=mtemp(end)+sum(mm(~tcMask));
                    massDist=mtemp./sum(mtemp);
                    
                    
                    massHist.inCGM.tcool(:,id+1)=xx;
                    massHist.inCGM.tcoolMass(:,id+1)=massDist;
                    
                    tCoolStruct.inCGM.meanTcMW(id+1)=sum(mm.*tc)/sum(mm);
                    tCoolStruct.inCGM.stdTcMW(id+1)=calc_standardDev(tmp,mm);
                    tCoolStruct.inCGM.medianTc(id+1)=10.^mus(3);
                    tCoolStruct.inCGM.quantTc(:,id+1)=10.^mus([1 2 4 5]);
                end
                % temp
                [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
                
                massHist.inCGM.temp(:,id+1)=xx;
                massHist.inCGM.tempMass(:,id+1)=massDist;
                
                tCoolStruct.inCGM.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inCGM.stdTempMW(id+1)=calc_standardDev(tmp,mm);
                tCoolStruct.inCGM.medianTemp(id+1)=10.^mus(3);
                tCoolStruct.inCGM.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
                
                %entropy
                [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
                
                massHist.inCGM.ent(:,id+1)=xx;
                massHist.inCGM.entMass(:,id+1)=massDist;
                
                tCoolStruct.inCGM.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inCGM.stdEntMW(id+1)=calc_standardDev(ent,mm);
                tCoolStruct.inCGM.medianEnt(id+1)=10.^mus(3);
                tCoolStruct.inCGM.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
                
                %number density
                [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
                
                massHist.inCGM.dens(:,id+1)=xx;
                massHist.inCGM.densMass(:,id+1)=massDist;
                
                tCoolStruct.inCGM.meanDensN(id+1)=mean(nDens);
                tCoolStruct.inCGM.stdDensN(id+1)=std(nDens);
                tCoolStruct.inCGM.medianDensN(id+1)=10.^mus(3);
                tCoolStruct.inCGM.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
                
                tCoolStruct.inCGM.gasMass(id+1)=sum(mm);
                tCoolStruct.inCGM.sfrMass(id+1)=sum(mass(distMask & sfMask));
                tCoolStruct.inCGM.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                
                tCoolStruct.inCGM.cellNum(id+1)=sum(mask);
                
                massHist.inCGM.gasMass(id+1)=sum(mm);
                
                
            end
            
            
            %% gas within 2* half stellar mass radius and edge
            distMask=gasDist>2.0.*rhalfStar ;
            
            mask=distMask & ~sfMask;
            
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);
                tmp=gas.Temperature(mask);
                ent=gas.Entropy(mask);
                nDens=nDensity(mask);
                
                
                %tcool
                tcMask=tc>0;
                if sum(tcMask)>0
                    [xx, massDist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus);
                    mtemp=massDist.*sum(mm(tcMask));
                    mtemp(end)=mtemp(end)+sum(mm(~tcMask));
                    massDist=mtemp./sum(mtemp);
                    
                    massHist.inSub.tcool(:,id+1)=xx;
                    massHist.inSub.tcoolMass(:,id+1)=massDist;
                    
                    tCoolStruct.inSub.meanTcMW(id+1)=sum(mm.*tc)/sum(mm);
                    tCoolStruct.inSub.stdTcMW(id+1)=calc_standardDev(tmp,mm);
                    tCoolStruct.inSub.medianTc(id+1)=10.^mus(3);
                    tCoolStruct.inSub.quantTc(:,id+1)=10.^mus([1 2 4 5]);
                end
                
                % temp
                [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
                
                massHist.inSub.temp(:,id+1)=xx;
                massHist.inSub.tempMass(:,id+1)=massDist;
                
                tCoolStruct.inSub.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inSub.stdTempMW(id+1)=calc_standardDev(tmp,mm);
                tCoolStruct.inSub.medianTemp(id+1)=10.^mus(3);
                tCoolStruct.inSub.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
                
                %entropy
                [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
                
                massHist.inSub.ent(:,id+1)=xx;
                massHist.inSub.entMass(:,id+1)=massDist;
                
                tCoolStruct.inSub.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inSub.stdEntMW(id+1)=calc_standardDev(ent,mm);
                tCoolStruct.inSub.medianEnt(id+1)=10.^mus(3);
                tCoolStruct.inSub.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
                
                %number density
                [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
                
                massHist.inSub.dens(:,id+1)=xx;
                massHist.inSub.densMass(:,id+1)=massDist;
                
                tCoolStruct.inSub.meanDensN(id+1)=mean(nDens);
                tCoolStruct.inSub.stdDensN(id+1)=std(nDens);
                tCoolStruct.inSub.medianDensN(id+1)=10.^mus(3);
                tCoolStruct.inSub.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
                
                
                tCoolStruct.inSub.gasMass(id+1)=sum(mm);
                tCoolStruct.inSub.sfrMass(id+1)=sum(mass(distMask & sfMask));
                tCoolStruct.inSub.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                
                tCoolStruct.inSub.cellNum(id+1)=sum(mask);
                
                massHist.inSub.gasMass(id+1)=sum(mm);
                
                
            end
            
            
                
            
            
            %% gas from  half gas mass radius and edge
            distMask=gasDist>rhalfGas ;
            
            mask=distMask & ~sfMask;
            
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);
                tmp=gas.Temperature(mask);
                ent=gas.Entropy(mask);
                nDens=nDensity(mask);
                
                
                %tcool
                tcMask=tc>0;
                if sum(tcMask)>0
                    [xx, massDist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus);
                    mtemp=massDist.*sum(mm(tcMask));
                    mtemp(end)=mtemp(end)+sum(mm(~tcMask));
                    massDist=mtemp./sum(mtemp);
                    
                    massHist.inOut.tcool(:,id+1)=xx;
                    massHist.inOut.tcoolMass(:,id+1)=massDist;
                    
                    tCoolStruct.inOut.meanTcMW(id+1)=sum(mm.*tc)/sum(mm);
                    tCoolStruct.inOut.stdTcMW(id+1)=calc_standardDev(tmp,mm);
                    tCoolStruct.inOut.medianTc(id+1)=10.^mus(3);
                    tCoolStruct.inOut.quantTc(:,id+1)=10.^mus([1 2 4 5]);
                end
                      
                % temp
                [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
                
                massHist.inOut.temp(:,id+1)=xx;
                massHist.inOut.tempMass(:,id+1)=massDist;
                
                tCoolStruct.inOut.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.inOut.stdTempMW(id+1)=calc_standardDev(tmp,mm);
                tCoolStruct.inOut.medianTemp(id+1)=10.^mus(3);
                tCoolStruct.inOut.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
                
                %entropy
                [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
                
                massHist.inOut.ent(:,id+1)=xx;
                massHist.inOut.entMass(:,id+1)=massDist;
                
                tCoolStruct.inOut.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.inOut.stdEntMW(id+1)=calc_standardDev(ent,mm);
                tCoolStruct.inOut.medianEnt(id+1)=10.^mus(3);
                tCoolStruct.inOut.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
                
                %number density
                [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
                
                massHist.inOut.dens(:,id+1)=xx;
                massHist.inOut.densMass(:,id+1)=massDist;
                
                tCoolStruct.inOut.meanDensN(id+1)=mean(nDens);
                tCoolStruct.inOut.stdDensN(id+1)=std(nDens);
                tCoolStruct.inOut.medianDensN(id+1)=10.^mus(3);
                tCoolStruct.inOut.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
                
                
                tCoolStruct.inOut.gasMass(id+1)=sum(mm);
                tCoolStruct.inOut.sfrMass(id+1)=sum(mass(distMask & sfMask));
                tCoolStruct.inOut.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                
                tCoolStruct.inOut.cellNum(id+1)=sum(mask);
                
                massHist.inOut.gasMass(id+1)=sum(mm);
                
            end
            
            
            
        end
        
        
        
        
        
    end
    
    
    %% save to mat file
    global DRACOFLAG
    if DRACOFLAG
        global DEFAULT_MATFILE_DIR
        global simDisplayName
        fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
        save([DEFAULT_MATFILE_DIR '/' fname],'tCoolStruct','-v7.3')
        
        fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
        
        fname2=sprintf('gasProperties_massHistograms_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
        save([DEFAULT_MATFILE_DIR '/' fname2],'massDist','-v7.3')
        
        fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname2]);
        
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
    
    
    
