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
global DEFAULT_MATFILE_DIR
global simDisplayName

if ~exist('readFlag','var')
    readFlag=true;
end

if readFlag
    fprintf(' *** Reading data *** \n');
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    readFlag=false;
    
    load([DEFAULT_MATFILE_DIR '/freeFallTime_profiles_snp' num2str(illUnits.snap) '_' simDisplayName '.mat']);
end


massThresh=10^9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
tCoolStruct.galMass=massAllGals;
tCoolStruct.galMask=galMask;
massHist.galMass=massAllGals;
massHist.galMask=galMask;

r200c=fofs.Group_R_Crit200(subsInfo.hostFof+1);

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;

len2=sum(galMask);
distLen=50;


%% initialize to zero

tcBin=[-1 0 1];
tctffBin=[-1 0 1];
tempBin=[5 6.5];
entBin=(-2:1:2);
densBin=(-4:1:-1);

tCoolStruct.tcBin=tcBin;
tCoolStruct.tctffBin=tctffBin;
tCoolStruct.tempBin=tempBin;
tCoolStruct.entBin=entBin;
tCoolStruct.densBin=densBin;


% tCoolStruct.inGal.meanTc=zeros(1,len);
tCoolStruct.inGal.meanTcMW=zeros(1,len);
tCoolStruct.inGal.medianTc=zeros(1,len);
tCoolStruct.inGal.modeTc=zeros(1,len);
tCoolStruct.inGal.quantTc=zeros(4,len);
tCoolStruct.inGal.massBinTc=zeros(length(tcBin)+1,len);
tCoolStruct.inGal.stdTcMW=zeros(1,len);

massHist.inGal.tc=zeros(distLen,len2);
massHist.inGal.tcMass=zeros(distLen,len2);

tCoolStruct.inGal.meanTcTffMW=zeros(1,len);
tCoolStruct.inGal.medianTcTff=zeros(1,len);
tCoolStruct.inGal.modeTcTff=zeros(1,len);
tCoolStruct.inGal.quantTcTff=zeros(4,len);
tCoolStruct.inGal.massBinTcTff=zeros(length(tctffBin)+1,len);
tCoolStruct.inGal.stdTcTffMW=zeros(1,len);

massHist.inGal.tcTff=zeros(distLen,len2);
massHist.inGal.tcTffMass=zeros(distLen,len2);

tCoolStruct.inGal.meanTempMW=zeros(1,len);
tCoolStruct.inGal.medianTemp=zeros(1,len);
tCoolStruct.inGal.modeTemp=zeros(1,len);
tCoolStruct.inGal.quantTemp=zeros(4,len);
tCoolStruct.inGal.massBinTemp=zeros(length(tempBin)+1,len);
tCoolStruct.inGal.stdTempMW=zeros(1,len);

massHist.inGal.temp=zeros(distLen,len2);
massHist.inGal.tempMass=zeros(distLen,len2);

tCoolStruct.inGal.meanEntMW=zeros(1,len);
tCoolStruct.inGal.medianEnt=zeros(1,len);
tCoolStruct.inGal.modeEnt=zeros(1,len);
tCoolStruct.inGal.quantEnt=zeros(4,len);
tCoolStruct.inGal.massBinEnt=zeros(length(entBin)+1,len);
tCoolStruct.inGal.stdEntMW=zeros(1,len);

massHist.inGal.ent=zeros(distLen,len2);
massHist.inGal.entMass=zeros(distLen,len2);

tCoolStruct.inGal.meanDensN=zeros(1,len);
tCoolStruct.inGal.medianDensN=zeros(1,len);
tCoolStruct.inGal.modeDensN=zeros(1,len);
tCoolStruct.inGal.quantDensN=zeros(4,len);
tCoolStruct.inGal.stdDensN=zeros(1,len);
tCoolStruct.inGal.massBinDensN=zeros(length(densBin)+1,len);

massHist.inGal.dens=zeros(distLen,len2);
massHist.inGal.densMass=zeros(distLen,len2);

tCoolStruct.inGal.gasMass=zeros(1,len);
tCoolStruct.inGal.sfrMass=zeros(1,len);
tCoolStruct.inGal.avgDensN=zeros(1,len);
tCoolStruct.inGal.cellNum=zeros(1,len);
tCoolStruct.inGal.sfr=zeros(1,len);

tCoolStruct.inCGM=tCoolStruct.inGal;
tCoolStruct.inSub=tCoolStruct.inGal;
tCoolStruct.inOut=tCoolStruct.inGal;

massHist.inCGM=massHist.inGal;
massHist.inSub=massHist.inGal;
massHist.inOut=massHist.inGal;

qus=[0.1 0.25 0.5 0.75 0.9];
tCoolStruct.qants=qus([1 2 4 5]);
cnt=0;
for id=0:len-1
    
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if galMask(id+1)
        
        
        cnt=cnt+1;
        % load gas from in sub halo
        
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas');
        %,{'Coordinates','Density','ElectronAbundance',...
        %    'InternalEnergy','Masses','GFM_CoolingRate','StarFormationRate'});%'EnergyDissipation','Machnumber','Velocities',...
        
        
        if gas.count==0
            continue
        end
        
        % get density and temperature
        nDensity=double(gas.Density.*illUnits.numberDensityFactor); %in cm^-3
        
        gas=illustris.utils.addTemperature( gas );
        %temp=illustris.utils.calcTemperature(gas.InternalEnergy,gas.ElectronAbundance); %in K
        
        
        mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
        
        
        if isfield(gas,'GFM_CoolingRate')
            tcool=illustris.utils.calcCoolingTime(gas.Density,gas.InternalEnergy,gas.GFM_CoolingRate); % cooling time in Gyr^-1
        else
            tcool=nan(size(mass));
        end
        
        
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
        
        % find free-fall time
        
        tff=recreate_TFF(gasDist,double(r200c(id+1)),...
            tffProfile.polyfit.a(id+1),...
            tffProfile.polyfit.b(id+1),...
            tffProfile.polyfit.c(id+1));
        tctff=tcool./tff;
        
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
        
        
        for kk=1:4
            
            switch kk
                case 1
                    % for gas within the gal (2*r_half - stellar)
                    distMask=gasDist<=2.0.*rhalfStar;
                    fld='inGal';
                    
                case 2
                    % for gas within the cgm (2*r_half  - r_half gas)
                    distMask=gasDist>2.0.*rhalfStar &  ...
                        gasDist<=rhalfGas;
                    fld='inCGM';
                    
                case 3
                    % gas from  half gas mass radius and edge
                    distMask=gasDist>rhalfGas ;
                    fld='inOut';
                    
                case 4
                    % gas within 2* half stellar mass radius and edge
                    distMask=gasDist>2.0.*rhalfStar ;
                    fld='inSub';
                    
            end
            
            % build mask 
            mask=distMask & ~sfMask;
            
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);
                tcff=tctff(mask);
                tmp=gas.Temperature(mask);
                ent=gas.Entropy(mask);
                nDens=nDensity(mask);
                
                
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                if sum(tcMask)>0
                    [xx, massDist, mus]=mk_mass_histogram(log10(tc(tcMask)),mm(tcMask),qus,distLen);
                    massDist(end)=massDist(end)+sum(mm(~tcMask));
                    [~,mxInd]=max(massDist);
                    
                    massHist.(fld).tc(:,cnt)=xx;
                    massHist.(fld).tcMass(:,cnt)=massDist;
                    
                    tCoolStruct.(fld).meanTcMW(id+1)=sum(mm(tcMask).*tc(tcMask))/sum(mm(tcMask));
                    tCoolStruct.(fld).stdTcMW(id+1)=calc_standardDev(tc(tcMask),mm(tcMask));
                    tCoolStruct.(fld).medianTc(id+1)=10.^mus(3);
                    tCoolStruct.(fld).quantTc(:,id+1)=10.^mus([1 2 4 5]);
                    tCoolStruct.(fld).modeTc(id+1)=xx(mxInd);
                                        
                    
                    ll=length(tcBin)+1;
                    for i=1:ll
                        if i==1
                            msk=xx<=tcBin(i);
                        elseif i==ll
                            msk=xx>tcBin(i-1);
                        else
                            msk=xx>tcBin(i-1) & xx<=tcBin(i);
                        end
                        
                        tCoolStruct.(fld).massBinTc(i,id+1)=sum(massDist(msk));
                    end
                    
                    % tc / tff
                    [xx, massDist, mus]=mk_mass_histogram(log10(tcff(tcMask)),mm(tcMask),qus,distLen);
                    massDist(end)=massDist(end)+sum(mm(~tcMask));
                                        [~,mxInd]=max(massDist);

                    massHist.(fld).tcTff(:,cnt)=xx;
                    massHist.(fld).tcTffMass(:,cnt)=massDist;
                    
                    tCoolStruct.(fld).meanTcTffMW(id+1)=sum(mm(tcMask).*tcff(tcMask))/sum(mm(tcMask));
                    tCoolStruct.(fld).stdTcTffMW(id+1)=calc_standardDev(tcff(tcMask),mm(tcMask));
                    tCoolStruct.(fld).medianTcTff(id+1)=10.^mus(3);
                    tCoolStruct.(fld).quantTcTff(:,id+1)=10.^mus([1 2 4 5]);
                    tCoolStruct.(fld).modeTcTff(id+1)=xx(mxInd);
                                        
                    ll=length(tctffBin)+1;
                    for i=1:ll
                        if i==1
                            msk=xx<=tctffBin(i);
                        elseif i==ll
                            msk=xx>tctffBin(i-1);
                        else
                            msk=xx>tctffBin(i-1) & xx<=tctffBin(i);
                        end
                        
                        tCoolStruct.(fld).massBinTcTff(i,id+1)=sum(massDist(msk));
                    end
                    
                    
                end
                
                % temp
                [xx, massDist, mus]=mk_mass_histogram(log10(tmp),mm,qus,distLen);
                                    [~,mxInd]=max(massDist);

                massHist.(fld).temp(:,cnt)=xx;
                massHist.(fld).tempMass(:,cnt)=massDist;
                
                tCoolStruct.(fld).meanTempMW(id+1)=sum(mm.*tmp)/sum(mm);
                tCoolStruct.(fld).stdTempMW(id+1)=calc_standardDev(tmp,mm);
                tCoolStruct.(fld).medianTemp(id+1)=10.^mus(3);
                tCoolStruct.(fld).quantTemp(:,id+1)=10.^mus([1 2 4 5]);
                             tCoolStruct.(fld).modeTemp(id+1)=xx(mxInd);   
                ll=length(tempBin)+1;
                for i=1:ll
                    if i==1
                        msk=xx<=tempBin(i);
                    elseif i==ll
                        msk=xx>tempBin(i-1);
                    else
                        msk=xx>tempBin(i-1) & xx<=tempBin(i);
                    end
                    
                    tCoolStruct.(fld).massBinTemp(i,id+1)=sum(massDist(msk));
                end
                
                %entropy
                [xx, massDist, mus]=mk_mass_histogram(log10(ent),mm,qus,distLen);
                                    [~,mxInd]=max(massDist);

                massHist.(fld).ent(:,cnt)=xx;
                massHist.(fld).entMass(:,cnt)=massDist;
                
                tCoolStruct.(fld).meanEntMW(id+1)=sum(mm.*ent)/sum(mm);
                tCoolStruct.(fld).stdEntMW(id+1)=calc_standardDev(ent,mm);
                tCoolStruct.(fld).medianEnt(id+1)=10.^mus(3);
                tCoolStruct.(fld).quantEnt(:,id+1)=10.^mus([1 2 4 5]);
                   tCoolStruct.(fld).modeEnt(id+1)=xx(mxInd);             
                   
                ll=length(entBin)+1;
                for i=1:ll
                    if i==1
                        msk=xx<=entBin(i);
                    elseif i==ll
                        msk=xx>entBin(i-1);
                    else
                        msk=xx>entBin(i-1) & xx<=entBin(i);
                    end
                    
                    tCoolStruct.(fld).massBinEnt(i,id+1)=sum(massDist(msk));
                end
                
                %number density
                [xx, massDist, mus]=mk_mass_histogram(log10(nDens),mm,qus,distLen);
                                    [~,mxInd]=max(massDist);

                massHist.(fld).dens(:,cnt)=xx;
                massHist.(fld).densMass(:,cnt)=massDist;
                
                tCoolStruct.(fld).meanDensN(id+1)=mean(nDens);
                tCoolStruct.(fld).stdDensN(id+1)=std(nDens);
                tCoolStruct.(fld).medianDensN(id+1)=10.^mus(3);
                tCoolStruct.(fld).quantDensN(:,id+1)=10.^mus([1 2 4 5]);
                tCoolStruct.(fld).modeDensN(id+1)=xx(mxInd);                                
                
                ll=length(densBin)+1;
                for i=1:ll
                    if i==1
                        msk=xx<=densBin(i);
                    elseif i==ll
                        msk=xx>densBin(i-1);
                    else
                        msk=xx>densBin(i-1) & xx<=densBin(i);
                    end
                    
                    tCoolStruct.(fld).massBinDensN(i,id+1)=sum(massDist(msk));
                end
                
                
                tCoolStruct.(fld).gasMass(id+1)=sum(mm); % only non-sf gas
                tCoolStruct.(fld).sfrMass(id+1)=sum(mass(distMask & sfMask));
                tCoolStruct.(fld).avgDensN(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
                tCoolStruct.(fld).sfr(id+1)=sum(gas.StarFormationRate(distMask & sfMask));
                                
                tCoolStruct.(fld).cellNum(id+1)=sum(mask);
                
                massHist.(fld).gasMass(cnt)=sum(mm);
                
                
            end
            
        end
        
    end
    
end


%% save to mat file
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fname=sprintf('gasProperties_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'tCoolStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
    fname2=sprintf('gasProperties_massHistograms_snp%s_%s',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname2],'massHist','-v7.3')
    
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



