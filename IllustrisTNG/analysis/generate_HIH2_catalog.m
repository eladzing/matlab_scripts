%% This Script generates HI and H_2 masses for eash subhalo - also breaking it down to CGM components

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
    
    HydroCat=illustris.utils.read_catalog('gas','snap',99,'folder','hydrogen');
    readFlag=false;
    
    
end

% assign a lower stellar mass threshold based on the TNG box (only TNG100)
massThresh=10^9; % threshold for *stellar* mass
fprintf('Setting stellar mass threshold to: %0.1e solar mass \n',massThresh);

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);

%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'gas');
hih2Struct.galMass=massAllGals;
hih2Struct.galMask=galMask;

%r200c=fofs.Group_R_Crit200(subsInfo.hostFof+1);



%% initialize to zero
fprintf('Initializing output... \n');

% set component
compNames=["Gal", "CGMin", "CGMout", "CGMall" "Sub"];
paramNames=["H2", "HI", "Gas", "Sfr"];
propNames="Mass";%, "avgDens"];

for fld=compNames
    hih2Struct.(fld).mask=int8(galMask);
    hih2Struct.(fld).Hmodel={'BR','GK','KMT'};
    
    for param=paramNames
        for prop=propNames
            pname=strcat(fld,param,prop);
            
            switch param
                case {"H2", "HI"}
                    hih2Struct.(fld).(pname)=zeros(3,len);
                otherwise
                    hih2Struct.(fld).(pname)=zeros(1,len);
            end
        end
    end
    
end


%% add generalities

%% run over subhalos
step=5;
stepNext=5;
len=double(subs.count);

len2=sum(galMask);

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
        fields={'Masses','Density','Coordinates','NeutralHydrogenAbundance','StarFormationRate'};
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas',fields);
        
        if gas.count==0
            continue
        end
        
        % get density and temperature
        %nDensity=double(gas.Density.*illUnits.numberDensityFactor); %in cm^-3
        
        %set mass
        mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
        
        % set H2 HI masses;
        subset = illustris.snapshot.getSnapOffsets(bp,snap,id,'Subhalo');% get the relevant particle indices
        %particle offests are:
        firstInd=subset.offsetType(1)+1;
        lastInd=firstInd+int64(subset.lenType(1))-1;
        
        gas.mH=HydroCat.MH(firstInd:lastInd)';
        gas.mH2(1,:)=HydroCat.MH2BR(firstInd:lastInd)';
        gas.mH2(2,:)=HydroCat.MH2GK(firstInd:lastInd)';
        gas.mH2(3,:)=HydroCat.MH2KMT(firstInd:lastInd)';
        gas.mHi(1,:)=gas.mH-HydroCat.MH2BR(firstInd:lastInd)';
        gas.mHi(2,:)=gas.mH-HydroCat.MH2GK(firstInd:lastInd)';
        gas.mHi(3,:)=gas.mH-HydroCat.MH2KMT(firstInd:lastInd)';
        
        
        % find distance from galaxy center
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
        gasDist=sqrt( sum(double(gas.newCoord).^2,1));
        
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
            
            % build mask -  within the ditance limitations
            mask=distMask ;
            
            if any(mask)
                % apply mask
                mm=mass(mask);
                %nDens=nDensity(mask);
                %mH=gas.mH(mask);
                
                for k=1:3
                    
                    param='HI';
                    hih2Struct.(fld).(strcat(fld,param,'mass'))(k,id+1)=sum(gas.mHi(k,mask));
                    param='H2';
                    hih2Struct.(fld).(strcat(fld,param,'mass'))(k,id+1)=sum(gas.mH2(k,mask));
                end
                
                
                
                hih2Struct.(fld).(strcat(fld, 'GasMass'))(id+1)=sum(mm); % All gas Mass
                hih2Struct.(fld).(strcat(fld, 'SfrMass'))(id+1)=sum(mass(distMask & sfMask));
                hih2Struct.(fld).([fld 'SFR'])(id+1)=sum(gas.StarFormationRate(distMask & sfMask));
                %hih2Struct.(fld).(strcat(fld, 'AvgDens'))(id+1)=sum(mm)/sum(mm./nDens); % in cm^-3
            end
            
        end
        
    end
end

%% save to catalog file
global DRACOFLAG
if DRACOFLAG
    
    
    
    %     for fld=compNames
    %
    %         outStruct=hih2Struct.(fld);
    %
    %         catName=sprintf('Subhalos_%s_PhysicalGasProperties',fld);
    %
    %         folder=['PhysicalGasProperties/' simDisplayName];
    %
    %         illustris.utils.write_catalog(outStruct,snap,'name',catName,...
    %             'path','default','folder','PhysicalGasProperties','v');
    %
    %
    %         %
    %         %         fname2=sprintf('gasProperties_massHistograms_snp%s_%s',num2str(snap),simDisplayName);
    %         %         save([DEFAULT_MATFILE_DIR '/' fname2],'massHist','-v7.3')
    %         %
    %         %         fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname2]);
    %         %
    %     end
    
    fname=sprintf('hih2Catalog_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'hih2Struct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
end
fprintf(' *** DONE!  *** \n');

