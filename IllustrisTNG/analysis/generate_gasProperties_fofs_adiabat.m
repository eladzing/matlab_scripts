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
%snap=99; %z=0
%bp=illustris.set_env(simName); %,'draco');

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
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
%% select galaxies

%massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
%galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals','hasGas');

galMask=(fofs.Group_M_Crit200.*illUnits.massUnit)>1e11;
%gasPropsFofs.galMass=massAllGals;
gasPropsFofs.galMask=galMask;

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(fofs.count);
cnt=0;

%distLen=200;


%% initialize to zero

tcBin=[0 1];
tempBin=[5 6.5];
entBin=(-2:1:2);
densBin=(-4:1:-1);

gasPropsFofs.tcBin=tcBin;
gasPropsFofs.tempBin=tempBin;
gasPropsFofs.entBin=entBin;
gasPropsFofs.densBin=densBin;


% gasPropsFofs.inR500.meanTc=zeros(1,len);
% gasPropsFofs.inR500.meanTcMW=zeros(1,len);
% gasPropsFofs.inR500.medianTc=zeros(1,len);
% gasPropsFofs.inR500.quantTc=zeros(4,len);
% gasPropsFofs.inR500.massBinTc=zeros(length(tcBin)+1,len);
% gasPropsFofs.inR500.stdTcMW=zeros(1,len);

%massHist.inR500.tcool=zeros(distLen,len);
%massHist.inR500.tcoolMass=zeros(distLen,len);

gasPropsFofs.inR500.meanTempMW=zeros(1,len);
gasPropsFofs.inR500.medianTemp=zeros(1,len);
gasPropsFofs.inR500.quantTemp=zeros(4,len);
gasPropsFofs.inR500.massBinTemp=zeros(length(tempBin)+1,len);
gasPropsFofs.inR500.stdTempMW=zeros(1,len);

%massHist.inR500.temp=zeros(distLen,len);
%massHist.inR500.tempMass=zeros(distLen,len);

gasPropsFofs.inR500.meanEntMW=zeros(1,len);
gasPropsFofs.inR500.medianEnt=zeros(1,len);
gasPropsFofs.inR500.quantEnt=zeros(4,len);
gasPropsFofs.inR500.massBinEnt=zeros(length(entBin)+1,len);
gasPropsFofs.inR500.stdEntMW=zeros(1,len);

%massHist.inR500.ent=zeros(distLen,len);
%massHist.inR500.entMass=zeros(distLen,len);

gasPropsFofs.inR500.meanDensN=zeros(1,len);
gasPropsFofs.inR500.medianDensN=zeros(1,len);
gasPropsFofs.inR500.quantDensN=zeros(4,len);
gasPropsFofs.inR500.stdDensN=zeros(1,len);
gasPropsFofs.inR500.massBinDens=zeros(length(densBin)+1,len);

%massHist.inR500.dens=zeros(distLen,len);
%massHist.inR500.densMass=zeros(distLen,len);

gasPropsFofs.inR500.gasMass=zeros(1,len);
%gasPropsFofs.inR500.sfrMass=zeros(1,len);
gasPropsFofs.inR500.avgDensN=zeros(1,len);
gasPropsFofs.inR500.cellNum=zeros(1,len);
%massHist.inR500.gasMass=zeros(1,len);


gasPropsFofs.inR200=gasPropsFofs.inR500;
gasPropsFofs.inTot=gasPropsFofs.inR500;


qus=[0.1 0.25 0.5 0.75 0.9];
gasPropsFofs.qants=qus([1 2 4 5]);

for id=0:len-1
    
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if galMask(id+1)
        
        
        % load gas form host of central galaxy
        idHost=id;   %subsInfo.hostFof(id+1);
        gas=illustris.snapshot.loadHalo(bp, snap, idHost, 'gas');
%         {'Coordinates','Density','ElectronAbundance',...
%             'InternalEnergy','Masses','GFM_CoolingRate','StarFormationRate'});%'EnergyDissipation','Machnumber','Velocities',...
%         
        
        if gas.count==0
            continue
        end
        
        nDensity=double(gas.Density.*illUnits.numberDensityFactor); %in cm^-3
        
        gas=illustris.utils.addTemperature( gas );
        %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
        
        
        mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
%         if isfield(gas,'GFM_CoolingRate')
%             tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr^-1
%         else
%             tcool=nan(size(mass));
%         end
        gas=illustris.utils.addEntropy( gas );
        
        
        % find distance from galaxy center
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,fofs.GroupPos(:,idHost+1));
        
        gasDist=sqrt( sum(double(gas.newCoord).^2,1));
        
        
        % find importat radii
        rmax=max(gasDist);
        % rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)); % stellar half mass radius
        % rhalfGas=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1)); % gas half mass radius
        
        
        % identify Star-forming gas
        
        %sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
        
        
        %% for gas within r500
        
        distMask=gasDist<=fofs.Group_R_Crit500(idHost+1);
        mask=distMask ;
        
        if any(mask)
            mm=mass(mask);
%             tc=tcool(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            nDens=nDensity(mask);
            
            
%             %tcool
%             % a small amount of cells have tc=0 set for positive cooling
%             % rates (heating). We exclude them from the histogram and add
%             % them to the last bin.
%             tcMask=tc>0;
%             if sum(tcMask)>0
%                 [xx, massDist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus);
%                 massDist(end)=massDist(end)+sum(mm(~tcMask));
%                 
%                 %massHist.inR500.tcool(:,id+1)=xx;
%                 %massHist.inR500.tcoolMass(:,id+1)=massDist;
%                 
%                 gasPropsFofs.inR500.meanTcMW(id+1)=sum(mm(tcMask).*tc(tcMask))/sum(mm(tcMask));
%                 gasPropsFofs.inR500.stdTcMW(id+1)=calc_standardDev(tc(tcMask),mm(tcMask));
%                 gasPropsFofs.inR500.medianTc(id+1)=10.^mus(3);
%                 gasPropsFofs.inR500.quantTc(:,id+1)=10.^mus([1 2 4 5]);
%                 
%                 ll=length(tcBin)+1;
%                 for i=1:ll
%                     if i==1
%                         msk=xx<=tcBin(i);
%                     elseif i==ll
%                         msk=xx>tcBin(i-1);
%                     else
%                         msk=xx>tcBin(i-1) & xx<=tcBin(i);
%                     end
%                     
%                     gasPropsFofs.inR500.massBinTc(i,id+1)=sum(massDist(msk));
%                 end
%                 
%             end
%             
            % temp
            [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
            
            %massHist.inR500.temp(:,id+1)=xx;
            %massHist.inR500.tempMass(:,id+1)=massDist;
            
            gasPropsFofs.inR500.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
            gasPropsFofs.inR500.stdTempMW(id+1)=calc_standardDev(tmp,mm);
            gasPropsFofs.inR500.medianTemp(id+1)=10.^mus(3);
            gasPropsFofs.inR500.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(tempBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=tempBin(i);
                elseif i==ll
                    msk=xx>tempBin(i-1);
                else
                    msk=xx>tempBin(i-1) & xx<=tempBin(i);
                end
                
                gasPropsFofs.inR500.massBinTemp(i,id+1)=sum(massDist(msk));
            end
            
            %entropy
            [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
            
            %massHist.inR500.ent(:,id+1)=xx;
            %massHist.inR500.entMass(:,id+1)=massDist;
            
            gasPropsFofs.inR500.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
            gasPropsFofs.inR500.stdEntMW(id+1)=calc_standardDev(ent,mm);
            gasPropsFofs.inR500.medianEnt(id+1)=10.^mus(3);
            gasPropsFofs.inR500.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
            
            
            ll=length(entBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=entBin(i);
                elseif i==ll
                    msk=xx>entBin(i-1);
                else
                    msk=xx>entBin(i-1) & xx<=entBin(i);
                end
                
                gasPropsFofs.inR500.massBinEnt(i,id+1)=sum(massDist(msk));
            end
            
            %number density
            [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
            
            %massHist.inR500.dens(:,id+1)=xx;
            %massHist.inR500.densMass(:,id+1)=massDist;
            
            gasPropsFofs.inR500.meanDensN(id+1)=mean(nDens);
            gasPropsFofs.inR500.stdDensN(id+1)=std(nDens);
            gasPropsFofs.inR500.medianDensN(id+1)=10.^mus(3);
            gasPropsFofs.inR500.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(densBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=densBin(i);
                elseif i==ll
                    msk=xx>densBin(i-1);
                else
                    msk=xx>densBin(i-1) & xx<=densBin(i);
                end
                
                gasPropsFofs.inR500.massBinDens(i,id+1)=sum(massDist(msk));
            end
            
            
            gasPropsFofs.inR500.gasMass(id+1)=sum(mm);
%             gasPropsFofs.inR500.sfrMass(id+1)=sum(mass(distMask & sfMask));
            gasPropsFofs.inR500.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
            
            gasPropsFofs.inR500.cellNum(id+1)=sum(mask);
            
            %massHist.inR500.gasMass(id+1)=sum(mm);
            
            
        end
        
        
        
        %% for gas within r200
        
        distMask=gasDist<=fofs.Group_R_Crit200(idHost+1);
        
        mask=distMask ;
        
        if any(mask)
            mm=mass(mask);
           
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            nDens=nDensity(mask);
            
            
            
            % temp
            [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
            
            %massHist.inR200.temp(:,id+1)=xx;
            %massHist.inR200.tempMass(:,id+1)=massDist;
            
            gasPropsFofs.inR200.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
            gasPropsFofs.inR200.stdTempMW(id+1)=calc_standardDev(tmp,mm);
            gasPropsFofs.inR200.medianTemp(id+1)=10.^mus(3);
            gasPropsFofs.inR200.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
            
            
            ll=length(tempBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=tempBin(i);
                elseif i==ll
                    msk=xx>tempBin(i-1);
                else
                    msk=xx>tempBin(i-1) & xx<=tempBin(i);
                end
                
                gasPropsFofs.inR200.massBinTemp(i,id+1)=sum(massDist(msk));
            end
            
            %entropy
            [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
            
            %massHist.inR200.ent(:,id+1)=xx;
            %massHist.inR200.entMass(:,id+1)=massDist;
            
            gasPropsFofs.inR200.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
            gasPropsFofs.inR200.stdEntMW(id+1)=calc_standardDev(ent,mm);
            gasPropsFofs.inR200.medianEnt(id+1)=10.^mus(3);
            gasPropsFofs.inR200.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(entBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=entBin(i);
                elseif i==ll
                    msk=xx>entBin(i-1);
                else
                    msk=xx>entBin(i-1) & xx<=entBin(i);
                end
                
                gasPropsFofs.inR200.massBinEnt(i,id+1)=sum(massDist(msk));
            end
            
            %number density
            [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
            
            %massHist.inR200.dens(:,id+1)=xx;
            %massHist.inR200.densMass(:,id+1)=massDist;
            
            gasPropsFofs.inR200.meanDensN(id+1)=mean(nDens);
            gasPropsFofs.inR200.stdDensN(id+1)=std(nDens);
            gasPropsFofs.inR200.medianDensN(id+1)=10.^mus(3);
            gasPropsFofs.inR200.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(densBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=densBin(i);
                elseif i==ll
                    msk=xx>densBin(i-1);
                else
                    msk=xx>densBin(i-1) & xx<=densBin(i);
                end
                
                gasPropsFofs.inR200.massBinDens(i,id+1)=sum(massDist(msk));
            end
            
            gasPropsFofs.inR200.gasMass(id+1)=sum(mm);
            gasPropsFofs.inR200.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
            
            gasPropsFofs.inR200.cellNum(id+1)=sum(mask);
            
            %massHist.inR200.gasMass(id+1)=sum(mm);
            
            
        end
        
        
        %% All the gas in the fof
        distMask=true(size(gasDist));
        
        mask=distMask;
        
        if any(mask)
            mm=mass(mask);
            tmp=gas.Temperature(mask);
            ent=gas.Entropy(mask);
            nDens=nDensity(mask);
            
           
            
            % temp
            [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus);
            
            %massHist.inTot.temp(:,id+1)=xx;
            %massHist.inTot.tempMass(:,id+1)=massDist;
            
            gasPropsFofs.inTot.meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
            gasPropsFofs.inTot.stdTempMW(id+1)=calc_standardDev(tmp,mm);
            gasPropsFofs.inTot.medianTemp(id+1)=10.^mus(3);
            gasPropsFofs.inTot.quantTemp(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(tempBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=tempBin(i);
                elseif i==ll
                    msk=xx>tempBin(i-1);
                else
                    msk=xx>tempBin(i-1) & xx<=tempBin(i);
                end
                
                gasPropsFofs.inTot.massBinTemp(i,id+1)=sum(massDist(msk));
            end
            
            %entropy
            [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus);
            
            %massHist.inTot.ent(:,id+1)=xx;
            %massHist.inTot.entMass(:,id+1)=massDist;
            
            gasPropsFofs.inTot.meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
            gasPropsFofs.inTot.stdEntMW(id+1)=calc_standardDev(ent,mm);
            gasPropsFofs.inTot.medianEnt(id+1)=10.^mus(3);
            gasPropsFofs.inTot.quantEnt(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(entBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=entBin(i);
                elseif i==ll
                    msk=xx>entBin(i-1);
                else
                    msk=xx>entBin(i-1) & xx<=entBin(i);
                end
                
                gasPropsFofs.inTot.massBinEnt(i,id+1)=sum(massDist(msk));
            end
            
            %number density
            [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus);
            
            %massHist.inTot.dens(:,id+1)=xx;
            %massHist.inTot.densMass(:,id+1)=massDist;
            
            gasPropsFofs.inTot.meanDensN(id+1)=mean(nDens);
            gasPropsFofs.inTot.stdDensN(id+1)=std(nDens);
            gasPropsFofs.inTot.medianDensN(id+1)=10.^mus(3);
            gasPropsFofs.inTot.quantDensN(:,id+1)=10.^mus([1 2 4 5]);
            
            ll=length(densBin)+1;
            for i=1:ll
                if i==1
                    msk=xx<=densBin(i);
                elseif i==ll
                    msk=xx>densBin(i-1);
                else
                    msk=xx>densBin(i-1) & xx<=densBin(i);
                end
                
                gasPropsFofs.inTot.massBinDens(i,id+1)=sum(massDist(msk));
            end
            
            gasPropsFofs.inTot.gasMass(id+1)=sum(mm);
            gasPropsFofs.inTot.avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
            
            gasPropsFofs.inTot.cellNum(id+1)=sum(mask);
            
            %massHist.inTot.gasMass(id+1)=sum(mm);
            
            
        end
        
        
    end
    
    
    
end








%% save to mat file
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    
    fname=sprintf('gasProperties_fofs_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'gasPropsFofs','-v7.3')
    
 
end
fprintf(' *** DONE!  *** \n');

%% save to postprocessed catalog

%global myPOSTPROCESSING
% gasPropsFofs.inR500.mask=floor(galMask);
% gasPropsFofs.inR200.mask=floor(galMask);
% gasPropsFofs.inTot.mask=floor(galMask);
% illustris.utils.write_catalog( gasPropsFofs.inR500,99,'name','gasProps_inR500','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( gasPropsFofs.inR200,99,'name','gasProps_inR200','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( gasPropsFofs.inTot,99,'name','gasProps_inTot','folderName','gasProperties','verbose');
%



