%% Study the differences between the gas in sf vs. quiescent centals
% go over all galaxies and extract cooling times for the gas and entropy as
% well
%
% several areas: within 2*r_h, all gas in the subhalo, others?

%% set framework

%global matFilePath
%global cosmoStruct
%global illUnits

%sim='300';
%snap=99; %z=0
%bp=illustris.set_env(simName,'draco');
illustris.utils.set_illUnits(snap);

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
illustris.utils.set_illUnits(snap)
global illUnits
%% selec galaxies

massAllGals= subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit; % stellar mass within 2*rhalf
massMask=massAllGals>massThresh;

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
%galMask= subs.SubhaloFlag & subsInfo.DM10perc & subsInfo.hasStars & massMask & subsInfo.hostHasVirial & subsInfo.hasGas ;%  & subs.SubhaloMassInRadType(1,:)>0;
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
%galMask=galMask & subsInfo.hasGas;

bhStruct.galMass=massAllGals;
bhStruct.galMask=galMask;

%% generate birds
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);


bhStruct.inGal.cumEngQM=zeros(1,len);
bhStruct.inGal.cumEngRM=zeros(1,len);

bhStruct.inGal.cumEngQM_maxBH=zeros(1,len);
bhStruct.inGal.cumEngRM_maxBH=zeros(1,len);

bhStruct.inGal.bhMassMax=zeros(1,len);
%bhStruct.inGal.bhMassCat=zeros(1,len);
bhStruct.inGal.bhMassSum=zeros(1,len);
bhStruct.inGal.bhCount=zeros(1,len);

bhStruct.inCGM=bhStruct.inGal;
bhStruct.inSub=bhStruct.inGal;
bhStruct.inOut=bhStruct.inGal;




for id=0:len-1
    
    
    perCent=floor((id+1)/len*100);
    
    if (perCent>=stepNext)     %       mod(perCent,5)==0 && perCent>=5)
        fprintf('%s %% of Galaxies done  \n',num2str(floor(perCent)))
        stepNext=stepNext+step;
    end
    
    if galMask(id+1)
        
        %cnt=cnt+1;
        
        
        
        % load bh data
        bhs=illustris.snapshot.loadSubhalo(bp, snap, id, 'bh'); %,{'Coordinates','Density','ElectronAbundance','InternalEnergy','Masses','GFM_CoolingRate'});
        
        
        % find distance from galaxy center
        if bhs.count>0
            
            
            bhs.newCoord = illustris.utils.centerObject(bhs.Coordinates,subs.SubhaloPos(:,id+1));
            bhsDist=sqrt( sum(bhs.newCoord.^2,1));
            
            
            % find importat radii
            %rmax=max(bhsDist);
            rhalfStar=subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1); % stellar half mass radius
            rhalfGas=subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1); % gas half mass radius
            %bhStruct.rmax(cnt)=rmax.*illUnits.lengthUnit;
            bhStruct.rhalfStar(id+1)=rhalfStar.*illUnits.lengthUnit;% in kpc
            bhStruct.rhalfGas(id+1)=rhalfGas.*illUnits.lengthUnit; % in kpc
            
            
            %     % identify Star-forming gas
            %     phaseDiagram_polygons
            %     sfMask=inpolygon(log10(dens),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
            %
            
            
            %% divide into differener gas components
            for kk=1:4
                
                switch kk
                    case 1
                        % for gas within the gal (2*r_half - stellar)
                        mask=bhsDist<=2.0.*rhalfStar;
                        fld='inGal';
                        
                    case 2
                        % for gas within the cgm (2*r_half  - r_half gas)
                        mask=bhsDist>2.0.*rhalfStar &  ...
                            bhsDist<=rhalfGas;
                        fld='inCGM';
                        
                    case 3
                        % gas from  half gas mass radius and edge
                        mask=bhsDist>rhalfGas ;
                        fld='inOut';
                        
                    case 4
                        % gas within 2* half stellar mass radius and edge
                        mask=bhsDist>2.0.*rhalfStar ;
                        fld='inSub';
                        
                end
                
                if any(mask)
                    
                    bhMass=bhs.BH_Mass(mask);
                    bhQM=bhs.BH_CumEgyInjection_QM(mask);
                    bhRM=bhs.BH_CumEgyInjection_RM(mask);
                    
                    [bhMax,maxInd]=max(bhs.BH_Mass(mask));
                    
                    bhStruct.(fld).bhMassMax(id+1)=bhMax.*illUnits.massUnit;
                    bhStruct.(fld).bhMassSum(id+1)=sum(bhMass).*illUnits.massUnit;
                    
                    %total
                    bhStruct.(fld).cumEngQM(id+1)=sum(bhQM).*illUnits.BHEnergyFactor;
                    bhStruct.(fld).cumEngRM(id+1)=sum(bhRM).*illUnits.BHEnergyFactor;
                    
                    %central
                    bhStruct.(fld).cumEngQM_maxBH(id+1)=bhQM(maxInd).*illUnits.BHEnergyFactor;
                    bhStruct.(fld).cumEngRM_maxBH(id+1)=bhRM(maxInd).*illUnits.BHEnergyFactor;
                    
                    bhStruct.(fld).bhCount(id+1)=sum(mask);
                    
                end
            end
            
            
            
            %
            %             %% for gas within the cgm (2*r_half  - r_half gas)
            %
            %             mask=bhsDist>2.0.*rhalfStar &  ...
            %                 bhsDist<=rhalfGas;
            %
            %
            %             if any(mask)
            %                 bhStruct.inCGM.cumEngQM(id+1)=sum(bhs.BH_CumEgyInjection_QM(mask)).*illUnits.BHEnergyFactor;
            %                 bhStruct.inCGM.cumEngRM(id+1)=sum(bhs.BH_CumEgyInjection_RM(mask)).*illUnits.BHEnergyFactor;
            %                 bhStruct.inCGM.bhMassMax(id+1)=max(bhs.BH_Mass(mask)).*illUnits.massUnit;
            %                 bhStruct.inCGM.bhMassSum(id+1)=sum(bhs.BH_Mass(mask)).*illUnits.massUnit;
            %                 bhStruct.inCGM.bhCount(id+1)=sum(mask);
            %
            %             end
            %
            %
            %
            %
            %
            %             %% gas within 2* half stellar mass radius and edge
            %             mask=bhsDist>2.0.*rhalfStar ;
            %
            %             if any(mask)
            %                 bhStruct.inSub.cumEngQM(id+1)=sum(bhs.BH_CumEgyInjection_QM(mask)).*illUnits.BHEnergyFactor;
            %                 bhStruct.inSub.cumEngRM(id+1)=sum(bhs.BH_CumEgyInjection_RM(mask)).*illUnits.BHEnergyFactor;
            %                 bhStruct.inSub.bhMassMax(id+1)=max(bhs.BH_Mass(mask)).*illUnits.massUnit;
            %                 bhStruct.inSub.bhMassSum(id+1)=sum(bhs.BH_Mass(mask)).*illUnits.massUnit;
            %                 bhStruct.inSub.bhCount(id+1)=sum(mask);
            %
            %             end
            
        end
    end
end


%% save to mat file
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    fname=sprintf('BH_energyInjection_snp%s_%s',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'bhStruct','-v7.3')
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
end

% bhStruct.inGal.mask=floor(galMask);
% bhStruct.inCGM.mask=floor(galMask);
% bhStruct.inSub.mask=floor(galMask);
% illustris.utils.write_catalog( bhStruct.inGal,99,'name','bh_inGal','folderName','bhProps');
% illustris.utils.write_catalog( bhStruct.inCGM,99,'name','bh_inCGM','folderName','bhProps');
% illustris.utils.write_catalog( bhStruct.inSub,99,'name','bh_inSub','folderName','bhProps');



fprintf(' *** DONE!  *** \n');



