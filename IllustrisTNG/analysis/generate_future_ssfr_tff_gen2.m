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
global cosmoStruct
global DEFAULT_MATFILE_DIR
global simDisplayName

if ~exist('readFlag','var')
    readFlag=true;
end

if readFlag
    fprintf(' *** Reading data *** \n');
    fofs=illustris.groupcat.loadHalos(bp,snap);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    
    load([DEFAULT_MATFILE_DIR '/freeFallTime_profiles_snp' num2str(snap) '_' simDisplayName '.mat'])
    
    readFlag=false;
end


massThresh=10^9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
%massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
%galMask= subs.SubhaloFlag & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;

galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
%galMask=galMask & subsInfo.hasGas;

sfrStruct.galMass=massAllGals;
sfrStruct.galMask=galMask;

sfrStruct.rhalfStar=zeros(size(galMask));
sfrStruct.rhalfGas=zeros(size(galMask));


%massHist.galMass=massAllGals;
%massHist.galMask=galMask;

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);
cnt=0;

%distLen=200;


%% initialize to zero

sfrStruct.inGal.sfr=zeros(1,len);

sfrStruct.inGal.fMs =zeros(1,len);
sfrStruct.inGal.fMs10 =zeros(1,len);


sfrStruct.inGal.gasMass=zeros(1,len);
sfrStruct.inGal.sfrMass=zeros(1,len);

sfrStruct.inGal.cellNum=zeros(1,len);
sfrStruct.inCGM=sfrStruct.inGal;
sfrStruct.inOut=sfrStruct.inGal;
sfrStruct.inSub=sfrStruct.inGal;

ts=1.0; %[0.1 0.5 1.0 1.33 1.6 2.0 5.0]; %in Gyr

sfeBase=0.01;

sfrStruct.ts=ts;
%
% ts100=100; % in Myr
% ts200=200; % in Myr
% ts500=500; % in Myr
% ts1000=1000; % in Myr - 1 Gyr
% ts100=100; % in Myr
% ts200=200; % in Myr
% ts500=500; % in Myr
% ts1000=1000; % in Myr - 1 Gyr
%
units;

for id=0:len-1
    
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if galMask(id+1)
        
        %% calculate free-fall profile
        
        
        
        %tff=sqrt(3*pi./(32*Units.G.*rhoMean))./Units.Gyr; % free-fall time in Gyr
        %tff=sqrt( 2*(rr.*Units.kpc).^3./(Units.G* mass(halo,rr,'kpc').*Units.Ms))./Units.Gyr;
        
        
        %cnt=cnt+1;
        % load gas from in sub halo
        
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas'...
            ,{'Coordinates','Density','ElectronAbundance',...
            'InternalEnergy','Masses','GFM_CoolingRate','StarFormationRate'});%'EnergyDissipation','Machnumber','Velocities',...
        
        
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
        
        % find distance from galaxy center
        
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
        gasDist=sqrt( sum(double(gas.newCoord).^2,1));
        r200c=fofs.Group_R_Crit200(subsInfo.hostFof(id+1)+1);
        
        %% find free-fall times at gas cell positions
        aCof=tffProfile.polyfit.a(id+1);
        bCof=tffProfile.polyfit.b(id+1);
        cCof=tffProfile.polyfit.c(id+1);
        
        tff=recreate_TFF(gasDist,r200c,aCof,bCof,cCof);
        
        tfFrac=(0:0.1:1.0).*ts;
        
        % find important radii
        rmax=max(gasDist);
        rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)); % stellar half mass radius
        rhalfGas=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1)); % gas half mass radius
        
        %tCoolStruct.rmax(id+)=rmax.*illUnits.lengthUnit;
        sfrStruct.rhalfStar(id+1)=rhalfStar.*illUnits.lengthUnit;% in kpc
        sfrStruct.rhalfGas(id+1)=rhalfGas.*illUnits.lengthUnit; % in kpc
        
        % identify Star-forming gas
        %polys=phaseDiagram_polygons('noplot');
        sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
        
        %% go over gas components
        for kk=1:4
            
            switch kk
                case 1
                    % for gas within the gal (2*r_half - stellar)
                    distMask=gasDist<=(2.0.*rhalfStar);
                    fld='inGal';
                    
                case 2
                    % for gas within the cgm (2*r_half  - r_half gas)
                    distMask=gasDist>(2.0.*rhalfStar) &  ...
                        gasDist<=rhalfGas;
                    fld='inCGM';
                    
                case 3
                    % gas from  half gas mass radius and edge
                    distMask=gasDist>rhalfGas ;
                    fld='inOut';
                    
                case 4
                    % gas within 2* half stellar mass radius and edge
                    distMask=gasDist>(2.0.*rhalfStar) ;
                    fld='inSub';
                    
            end
            
            sfr=sum(gas.StarFormationRate(distMask & sfMask));
            sfrStruct.(fld).sfr(id+1)=sfr;
            
            mask=distMask & ~sfMask;
            
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);
                tf=tff(mask);
                
%                 sfr=sum(gas.StarFormationRate(distMask & sfMask));
%                 sfrStruct.(fld).sfr(id+1)=sfr;
%                 
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                
                tcMask= tc<tf & tc>0;
                tcMask10= tc<(10.*tf) & tc>0;
                
                accumMass=0;
                if sum(tcMask)>0
                    %for k=1:length(ts)
                    for ii=2:length(tfFrac)
                        tfMask=tf>=tfFrac(ii-1) & tf<tfFrac(ii);
                        
                        q=(ts/tfFrac(ii));
                        sfe=1-(1-sfeBase).^q;
                        
                        accumMass=accumMass+sfe.*sum(mm(tcMask & tfMask));
                    end
                    sfrStruct.(fld).fMs(id+1)= accumMass ; %sum(mm(tcMask & tf<(ts(k))));
                    %end
                end
                
                accumMass=0;
                if sum(tcMask10)>0
                    % for k=1:length(ts)
                    for ii=2:length(tfFrac)
                        tfMask=tf>=tfFrac(ii-1) & tf<tfFrac(ii);
                        q=(ts/tfFrac(ii));
                        sfe=1-(1-sfeBase).^q;
                        accumMass=accumMass+sfe.*sum(mm(tcMask10 & tfMask));
                    end
                    
                    sfrStruct.(fld).fMs10(id+1)= accumMass; %sum(mm(tcMask10 & tf<(ts(k))));
                    % end
                end
                
                sfrStruct.(fld).gasMass(id+1)=sum(mm); % only non-sf gas
                
                
            end
            sfrStruct.(fld).sfrMass(id+1)=sum(mass(distMask & sfMask));
            sfrStruct.(fld).cellNum(id+1)=sum(mask);
            
        end
        
        
    end
end


%% save to mat file
global DRACOFLAG
if DRACOFLAG
    %         global DEFAULT_MATFILE_DIR
    %         global simDisplayName
    %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fname=sprintf('futureSFR_tff_gen2_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'sfrStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    %
    %         fname2=sprintf('gasProperties_massHistograms_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    %         save([DEFAULT_MATFILE_DIR '/' fname2],'massDist','-v7.3')
    %
    %         fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname2]);
    %
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



