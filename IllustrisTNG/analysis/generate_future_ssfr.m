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
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs,snap);
units; % load general unit structure in cgs.
%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
%massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
%galMask= subs.SubhaloFlag & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
%galMask=galMask & subsInfo.hasGas;
sfrStruct.galMass=massAllGals;
sfrStruct.galMask=galMask;
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

sfrStruct.inGal.fMs100MyrSFR =zeros(1,len);
sfrStruct.inGal.fMs100MyrCool=zeros(1,len);

sfrStruct.inGal.fMs200MyrSFR =zeros(1,len);
sfrStruct.inGal.fMs200MyrCool=zeros(1,len);

sfrStruct.inGal.fMs500MyrSFR =zeros(1,len);
sfrStruct.inGal.fMs500MyrCool= zeros(1,len);

sfrStruct.inGal.fMs1GyrSFR = zeros(1,len);
sfrStruct.inGal.fMs1GyrCool=zeros(1,len);


sfrStruct.inGal.gasMass=zeros(1,len);
sfrStruct.inGal.sfrMass=zeros(1,len);

sfrStruct.inGal.cellNum=zeros(1,len);

sfrStruct.inCGM=sfrStruct.inGal;
sfrStruct.inOut=sfrStruct.inGal;
sfrStruct.inSub=sfrStruct.inGal;


ts100=100; % in Myr
ts200=200; % in Myr
ts500=500; % in Myr
ts1000=1000; % in Myr - 1 Gyr


    for id=0:len-1
        
        
        perCent=floor((id+1)/len*100);
        
        if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
            fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
            stepNext=stepNext+step;
        end
        
        if galMask(id+1)
            
            
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
%                 tmp=gas.Temperature(mask);
%                 ent=gas.Entropy(mask);
%                 nDens=nDensity(mask);
%                                 
                sfr=sum(gas.StarFormationRate(distMask & sfMask));
                sfrStruct.inGal.sfr(id+1)=sfr;
                
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                
                sfrStruct.inGal.fMs100MyrSFR(id+1) = sfr.*ts100.*1e6;
                sfrStruct.inGal.fMs200MyrSFR(id+1) = sfr.*ts200.*1e6;
                sfrStruct.inGal.fMs500MyrSFR(id+1) = sfr.*ts500.*1e6;
                sfrStruct.inGal.fMs1GyrSFR(id+1)   = sfr.*ts1000.*1e6;
                
                if sum(tcMask)>0
                    
                    sfrStruct.inGal.fMs100MyrCool(id+1)= sum(mm(tcMask & tc<(ts100./1e3)));
                    
                    sfrStruct.inGal.fMs200MyrCool(id+1)= sum(mm(tcMask & tc<(ts200./1e3)));
                    
                    sfrStruct.inGal.fMs500MyrCool(id+1)= sum(mm(tcMask & tc<(ts500./1e3)));
                    
                    sfrStruct.inGal.fMs1GyrCool(id+1)= sum(mm(tcMask & tc<(ts1000./1e3)));
                    
                end
                    
                sfrStruct.inGal.gasMass(id+1)=sum(mm); % only non-sf gas
                sfrStruct.inGal.sfrMass(id+1)=sum(mass(distMask & sfMask));
                
                sfrStruct.inGal.cellNum(id+1)=sum(mask);
                
            end
            
            
            
            %% for gas within the cgm (2*r_half  - r_half gas)
            
            distMask=gasDist>2.0.*rhalfStar &  ...
                gasDist<=rhalfGas;
                        
            mask=distMask & ~sfMask;
                        
            if any(mask)
                mm=mass(mask);
                tc=tcool(mask);

                sfr=sum(gas.StarFormationRate(distMask & sfMask));
                sfrStruct.inCGM.sfr(id+1)=sfr;
                
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                
                sfrStruct.inCGM.fMs100MyrSFR(id+1) = sfr.*ts100.*1e6;
                sfrStruct.inCGM.fMs200MyrSFR(id+1) = sfr.*ts200.*1e6;
                sfrStruct.inCGM.fMs500MyrSFR(id+1) = sfr.*ts500.*1e6;
                sfrStruct.inCGM.fMs1GyrSFR(id+1)   = sfr.*ts1000.*1e6;
                
                if sum(tcMask)>0
                    
                    sfrStruct.inCGM.fMs100MyrCool(id+1)= sum(mm(tcMask & tc<(ts100./1e3)));
                    
                    sfrStruct.inCGM.fMs200MyrCool(id+1)= sum(mm(tcMask & tc<(ts200./1e3)));
                    
                    sfrStruct.inCGM.fMs500MyrCool(id+1)= sum(mm(tcMask & tc<(ts500./1e3)));
                    
                    sfrStruct.inCGM.fMs1GyrCool(id+1)= sum(mm(tcMask & tc<(ts1000./1e3)));
                    
                end
                
                
                                
                sfrStruct.inCGM.gasMass(id+1)=sum(mm); % only non-sf gas
                sfrStruct.inCGM.sfrMass(id+1)=sum(mass(distMask & sfMask));
                
                sfrStruct.inCGM.cellNum(id+1)=sum(mask);
                
            end
            
            
            
            %% gas within 2* half stellar mass radius and edge
            distMask=gasDist>2.0.*rhalfStar ;
            
            mask=distMask & ~sfMask;
            
            if any(mask)
                  mm=mass(mask);
                tc=tcool(mask);

                sfr=sum(gas.StarFormationRate(distMask & sfMask));
                sfrStruct.inSub.sfr(id+1)=sfr;
                
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                
                sfrStruct.inSub.fMs100MyrSFR(id+1) = sfr.*ts100.*1e6;
                sfrStruct.inSub.fMs200MyrSFR(id+1) = sfr.*ts200.*1e6;
                sfrStruct.inSub.fMs500MyrSFR(id+1) = sfr.*ts500.*1e6;
                sfrStruct.inSub.fMs1GyrSFR(id+1)   = sfr.*ts1000.*1e6;
                
                if sum(tcMask)>0
                    
                    sfrStruct.inSub.fMs100MyrCool(id+1)= sum(mm(tcMask & tc<(ts100./1e3)));
                    
                    sfrStruct.inSub.fMs200MyrCool(id+1)= sum(mm(tcMask & tc<(ts200./1e3)));
                    
                    sfrStruct.inSub.fMs500MyrCool(id+1)= sum(mm(tcMask & tc<(ts500./1e3)));
                    
                    sfrStruct.inSub.fMs1GyrCool(id+1)= sum(mm(tcMask & tc<(ts1000./1e3)));
                    
                end
                                
                sfrStruct.inSub.gasMass(id+1)=sum(mm); % only non-sf gas
                sfrStruct.inSub.sfrMass(id+1)=sum(mass(distMask & sfMask));
                
                sfrStruct.inSub.cellNum(id+1)=sum(mask);
            end
            
            
                
            
            
            %% gas from  half gas mass radius and edge
            distMask=gasDist>rhalfGas ;
            
            mask=distMask & ~sfMask;
            
            if any(mask)
                  mm=mass(mask);
                tc=tcool(mask);

                sfr=sum(gas.StarFormationRate(distMask & sfMask));
                sfrStruct.inOut.sfr(id+1)=sfr;
                
                %tcool
                % a small amount of cells have tc=0 set for positive cooling
                % rates (heating). We exclude them from the histogram and add
                % them to the last bin.
                tcMask=tc>0;
                
                sfrStruct.inOut.fMs100MyrSFR(id+1) = sfr.*ts100.*1e6;
                sfrStruct.inOut.fMs200MyrSFR(id+1) = sfr.*ts200.*1e6;
                sfrStruct.inOut.fMs500MyrSFR(id+1) = sfr.*ts500.*1e6;
                sfrStruct.inOut.fMs1GyrSFR(id+1)   = sfr.*ts1000.*1e6;
                
                if sum(tcMask)>0
                    
                    sfrStruct.inOut.fMs100MyrCool(id+1)= sum(mm(tcMask & tc<(ts100./1e3)));
                    
                    sfrStruct.inOut.fMs200MyrCool(id+1)= sum(mm(tcMask & tc<(ts200./1e3)));
                    
                    sfrStruct.inOut.fMs500MyrCool(id+1)= sum(mm(tcMask & tc<(ts500./1e3)));
                    
                    sfrStruct.inOut.fMs1GyrCool(id+1)= sum(mm(tcMask & tc<(ts1000./1e3)));
                    
                end
                                
                sfrStruct.inOut.gasMass(id+1)=sum(mm); % only non-sf gas
                sfrStruct.inOut.sfrMass(id+1)=sum(mass(distMask & sfMask));
                
                sfrStruct.inOut.cellNum(id+1)=sum(mask);
            end
            
            
            
        end
        
        
        
        
        
    end
    
    
    %% save to mat file
    global DRACOFLAG
    if DRACOFLAG
        global DEFAULT_MATFILE_DIR
        global simDisplayName
        %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
        fname=sprintf('futureSFR_snp%s_%s.mat',num2str(snap),simDisplayName);
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
    
    
    
