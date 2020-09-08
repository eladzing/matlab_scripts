%% This Script generates interesting gas properties for the subfind catalogs
% and saves them as catalogs for TNG


%% set framework
% the simulation, snapshot, and environment should be set prior to running
% the script

illustris.utils.set_illUnits(snap);

global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName

if ~exist('readFlag','var')
    readFlag=true;
end

% read FOF and SUBFIND data, as well as free-fall time profiles
if readFlag
    fprintf(' *** Reading data *** \n');
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
    
    load([DEFAULT_MATFILE_DIR '/freeFallTime_profiles_snp' num2str(illUnits.snap) '_' simDisplayName '.mat']);
end

% assign a lower stellar mass threshold based on the TNG box.
if contains(simDisplayName,'100') || contains(simDisplayName,'300')
    massThresh=10^9; % threshold for *stellar* mass
elseif contains(simDisplayName,'50')
    massThresh=10^7; % threshold for *stellar* mass
else
    error('mass threshold not set  - could not identify simulation: %s \n',simDisplayName);
end
fprintf('Setting stellar mass threshold to: %0.1e solar mass \n',massThresh);



%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.

%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'gas');
PropStruct.galMass=massAllGals;
PropStruct.galMask=galMask;

r200c=fofs.Group_R_Crit200(subsInfo.hostFof+1);

%% generate values
%fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);

len2=sum(galMask);

% size of mass histogram
distLen=100;


%% initialize to zero
fprintf('Initializing output... \n');

% set component
compNames=["Gal", "CGMin", "CGMout", "CGMall" "Sub"];
paramNames=["Tcool", "TcTff", "Temp", "Entropy", "Density"];
propNames=["MeanMW", "StdDevMW", "MassMedian", "MassQuantiles"];

% These values are not needed for the catalogs
% tcBin=[-1 0 1];
% tctffBin=[-1 0 1];
% tempBin=[5 6.5];
% entBin=(-2:1:2);
% densBin=(-4:1:-1);

% PropStruct.tcBin=tcBin;
% PropStruct.tctffBin=tctffBin;
% PropStruct.tempBin=tempBin;
% PropStruct.entBin=entBin;
% PropStruct.densBin=densBin;


%% initialize output

for fld=compNames
    PropStruct.(fld).mask=int8(galMask);
    for param=paramNames
        for prop=propNames
            pname=strcat(fld,param,prop);
            switch prop
                case 'MassQuantiles'
                    PropStruct.(fld).(pname)=zeros(4,len);
                otherwise
                    PropStruct.(fld).(pname)=zeros(1,len);
            end
        end
    end
    
    PropStruct.(fld).(strcat(fld,'GasMass'))=zeros(1,len);
    PropStruct.(fld).(strcat(fld,'SfrMass'))=zeros(1,len);
    PropStruct.(fld).(strcat(fld,'AvgDens'))=zeros(1,len);
end


qus=[0.1 0.25 0.5 0.75 0.9];
PropStruct.qants=qus([1 2 4 5]);

%% run over subind objects
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

cnt=0;
for id=0:len-1
    
    % for following progress
    perCent=floor((cnt+1)/len2*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    % throw away unneeded objects
    if galMask(id+1)
        
        
        cnt=cnt+1;
        % load gas from in sub halo
        
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas');
        
        if gas.count==0
            continue
        end
        
        % get density and temperature
        nDensity=double(gas.Density.*illUnits.numberDensityFactor); %in cm^-3
        
        gas=illustris.utils.addTemperature( gas );
        %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
        
        %set mass
        mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
        
        % calculate cooling time
        if isfield(gas,'GFM_CoolingRate')
            tcool=double(illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate)); % cooling time in Gyr^-1
        else
            tcool=nan(size(mass));
        end
        
        % add entropy
        gas=illustris.utils.addEntropy( gas );
        
        
        % find distance from galaxy center
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
        gasDist=sqrt( sum(double(gas.newCoord).^2,1));
        
        % find free-fall time
        tff=recreate_TFF(gasDist,double(r200c(id+1)),...
            tffProfile.polyfit.a(id+1),...
            tffProfile.polyfit.b(id+1),...
            tffProfile.polyfit.c(id+1));
        tctff=tcool./tff;
        
        % find important radii
        rmax=max(gasDist);
        rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)); % stellar half mass radius
        rhalfGas=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1)); % gas half mass radius
        
        % identify Star-forming gas
        sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
        
        %% go over components
        for fld=compNames
            
            switch fld
                case 'Gal'
                    % for gas within the gal (2*r_half - stellar)
                    distMask=gasDist<=2.0.*rhalfStar;
                    
                case 'CGMin'
                    % for gas within the cgm (2*r_half  - r_half gas)
                    distMask=gasDist>2.0.*rhalfStar &  ...
                        gasDist<=rhalfGas;
                    
                case 'CGMout'
                    % gas from  half gas mass radius and edge
                    distMask=gasDist>rhalfGas ;
                    
                case 'CGMall'
                    % gas within 2* half stellar mass radius and edge
                    distMask=gasDist>2.0.*rhalfStar ;
                    
                case 'Sub'
                    % all gas in the subfind object
                    distMask=true(size(gasDist));
            end
            
            % build mask - non-sfr and within the ditance limitations
            mask=distMask & ~sfMask;
            
            if any(mask)
                % apply mask
                mm=mass(mask);
                tc=tcool(mask);
                tcff=tctff(mask);
                tmp=gas.Temperature(mask);
                ent=gas.Entropy(mask);
                nDens=nDensity(mask);
                
                
                %% tcool and tc/tff
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                if sum(tcMask)>0
                    [~, ~, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus,distLen);
                    %massDist(end)=massDist(end)+sum(mm(~tcMask));
                    %[~,mxInd]=max(massDist);
                    param='Tcool';
                    PropStruct.(fld).(strcat(fld,param,'MeanMW'))(id+1)=sum(mm(tcMask).*tc(tcMask))/sum(mm(tcMask));
                    PropStruct.(fld).(strcat(fld,param,'StdDevMW'))(id+1)=calc_standardDev(tc(tcMask),mm(tcMask));
                    PropStruct.(fld).(strcat(fld,param,'MassMedian'))(id+1)=10.^mus(3);
                    PropStruct.(fld).(strcat(fld,param,'MassQuantiles'))(:,id+1)=10.^mus([1 2 4 5]);
                    %PropStruct.(fld).modeTc(id+1)=xx(mxInd);
                    
                    %                     ll=length(tcBin)+1;
                    %                     for i=1:ll
                    %                         if i==1
                    %                             msk=xx<=tcBin(i);
                    %                         elseif i==ll
                    %                             msk=xx>tcBin(i-1);
                    %                         else
                    %                             msk=xx>tcBin(i-1) & xx<=tcBin(i);
                    %                         end
                    %
                    %                         PropStruct.(fld).massBinTc(i,id+1)=sum(massDist(msk));
                    %                     end
                    
                    % tc / tff
                    [~, ~, mus]=mk_mass_histogram(log10(tcff(tcMask)),mm(tcMask),qus,distLen);
                    %massDist(end)=massDist(end)+sum(mm(~tcMask));
                    %[~,mxInd]=max(massDist);
                    
                    param='TcTff';
                    PropStruct.(fld).(strcat(fld,param,'MeanMW'))(id+1)=sum(mm(tcMask).*tcff(tcMask))/sum(mm(tcMask));
                    PropStruct.(fld).(strcat(fld,param,'StdDevMW'))(id+1)=calc_standardDev(tcff(tcMask),mm(tcMask));
                    PropStruct.(fld).(strcat(fld,param,'MassMedian'))(id+1)=10.^mus(3);
                    PropStruct.(fld).(strcat(fld,param,'MassQuantiles'))(:,id+1)=10.^mus([1 2 4 5]);
                    %PropStruct.(fld).modeTcTff(id+1)=xx(mxInd);
                    
                    %                     ll=length(tctffBin)+1;
                    %                     for i=1:ll
                    %                         if i==1
                    %                             msk=xx<=tctffBin(i);
                    %                         elseif i==ll
                    %                             msk=xx>tctffBin(i-1);
                    %                         else
                    %                             msk=xx>tctffBin(i-1) & xx<=tctffBin(i);
                    %                         end
                    %
                    %                         PropStruct.(fld).massBinTcTff(i,id+1)=sum(massDist(msk));
                    %                     end
                    
                    
                end
                
                %% temperature
                [~, ~, mus]=mk_mass_histogram(log10(tmp),mm,qus,distLen);
                %[~,mxInd]=max(massDist);
                param='Temp';
                PropStruct.(fld).(strcat(fld,param,'MeanMW'))(id+1)=sum(mm.*tmp)/sum(mm);
                PropStruct.(fld).(strcat(fld,param,'StdDevMW'))(id+1)=calc_standardDev(tmp,mm);
                PropStruct.(fld).(strcat(fld,param,'MassMedian'))(id+1)=10.^mus(3);
                PropStruct.(fld).(strcat(fld,param,'MassQuantiles'))(:,id+1)=10.^mus([1 2 4 5]);
                %PropStruct.(fld).modeTemp(id+1)=xx(mxInd);
                %
                %                 ll=length(tempBin)+1;
                %                 for i=1:ll
                %                     if i==1
                %                         msk=xx<=tempBin(i);
                %                     elseif i==ll
                %                         msk=xx>tempBin(i-1);
                %                     else
                %                         msk=xx>tempBin(i-1) & xx<=tempBin(i);
                %                     end
                %
                %                     PropStruct.(fld).massBinTemp(i,id+1)=sum(massDist(msk));
                %                 end
                
                %% entropy
                [~, ~, mus]=mk_mass_histogram(log10(ent),mm,qus,distLen);
                %[~,mxInd]=max(massDist);
                param='Entropy';
                PropStruct.(fld).(strcat(fld,param,'MeanMW'))(id+1)=sum(mm.*ent)/sum(mm);
                PropStruct.(fld).(strcat(fld,param,'StdDevMW'))(id+1)=calc_standardDev(ent,mm);
                PropStruct.(fld).(strcat(fld,param,'MassMedian'))(id+1)=10.^mus(3);
                PropStruct.(fld).(strcat(fld,param,'MassQuantiles'))(:,id+1)=10.^mus([1 2 4 5]);
                %PropStruct.(fld).modeEnt(id+1)=xx(mxInd);
                
                %                 ll=length(entBin)+1;
                %                 for i=1:ll
                %                     if i==1
                %                         msk=xx<=entBin(i);
                %                     elseif i==ll
                %                         msk=xx>entBin(i-1);
                %                     else
                %                         msk=xx>entBin(i-1) & xx<=entBin(i);
                %                     end
                %
                %                     PropStruct.(fld).massBinEnt(i,id+1)=sum(massDist(msk));
                %                 end
                
                %% number density
                [~, ~, mus]=mk_mass_histogram(log10(nDens),mm,qus,distLen);
                %[~,mxInd]=max(massDist);
                param='Density';
                PropStruct.(fld).(strcat(fld,param,'MeanMW'))(id+1)=mean(nDens);
                PropStruct.(fld).(strcat(fld,param,'StdDevMW'))(id+1)=std(nDens);
                PropStruct.(fld).(strcat(fld,param,'MassMedian'))(id+1)=10.^mus(3);
                PropStruct.(fld).(strcat(fld,param,'MassQuantiles'))(:,id+1)=10.^mus([1 2 4 5]);
                %PropStruct.(fld).modeDensN(id+1)=xx(mxInd);
                
                %                 ll=length(densBin)+1;
                %                 for i=1:ll
                %                     if i==1
                %                         msk=xx<=densBin(i);
                %                     elseif i==ll
                %                         msk=xx>densBin(i-1);
                %                     else
                %                         msk=xx>densBin(i-1) & xx<=densBin(i);
                %                     end
                %
                %                     PropStruct.(fld).massBinDensN(i,id+1)=sum(massDist(msk));
                %                 end
                %
                
                PropStruct.(fld).(strcat(fld, 'GasMass'))(id+1)=sum(mm); % only non-sf gas
                PropStruct.(fld).(strcat(fld, 'SfrMass'))(id+1)=sum(mass(distMask & sfMask));
                PropStruct.(fld).(strcat(fld, 'AvgDens'))(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                %PropStruct.(fld).([fld 'Sfr'])(id+1)=sum(gas.StarFormationRate(distMask & sfMask));
                %
                %                 PropStruct.(fld).cellNum(id+1)=sum(mask);
                %
            end
            
        end
        
    end
    
end


%% save to catalog file
global DRACOFLAG
if DRACOFLAG
    
    
    
    for fld=compNames
        
        outStruct=PropStruct.(fld);
        
        catName=sprintf('Subhalos_%s_PhysicalGasProperties',fld);
        
        folder=['PhysicalGasProperties/' simDisplayName];
        
        illustris.utils.write_catalog(outStruct,snap,'name',catName,...
            'path','default','folder','PhysicalGasProperties','v');
        
        
        %
        %         fname2=sprintf('gasProperties_massHistograms_snp%s_%s',num2str(snap),simDisplayName);
        %         save([DEFAULT_MATFILE_DIR '/' fname2],'massHist','-v7.3')
        %
        %         fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname2]);
        %
    end
    
    fname=sprintf('gasProperties_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'PropStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
end
fprintf(' *** DONE!  *** \n');

%% save to postprocessed catalog

%global myPOSTPROCESSING
% PropStruct.inGal.mask=floor(galMask);
% PropStruct.inCGM.mask=floor(galMask);
% PropStruct.inSub.mask=floor(galMask);
% illustris.utils.write_catalog( PropStruct.inGal,99,'name','gasProps_inGal','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( PropStruct.inCGM,99,'name','gasProps_inCGM','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( PropStruct.inSub,99,'name','gasProps_inSub','folderName','gasProperties','verbose');
%



