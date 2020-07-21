%% find the average cell size of gas cells
%  also with the division into SF & non-star forming cells

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
global DRACOFLAG

if ~exist('readFlag','var')
    readFlag=true;
end

if readFlag
    fprintf(' *** Reading data *** \n');
    
    if DRACOFLAG
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
        
    else
        loadFofSubTNG100
    end
    
    readFlag=false;
    
end


massThresh=10^9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

% find back-splash
if snap==99
spbMask = mk_splashback_mask('time',5,'both',0.1);

else
    spbMask = false(size(tCoolStruct.galMask));
end


% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
galMaskCentrals=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals'); 
galMaskSats=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'sats');
galMaskBksp=spbMask;


cellSizeStruct.galMass=massAllGals;
cellSizeStruct.galMask=galMask;
cellSizeStruct.galMaskCentrals=galMaskCentrals;
cellSizeStruct.galMaskSats=galMaskSats;
cellSizeStruct.galMaskBksp=spbMask;

%r200c=fofs.Group_R_Crit200(subsInfo.hostFof+1);

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);

len2=sum(galMask);
distLen=50;


%% initialize to zero


cellSizeStruct.inGal.All.clMean=zeros(1,len);
cellSizeStruct.inGal.All.clMean2=zeros(1,len);
cellSizeStruct.inGal.All.clMedian=zeros(1,len);
%cellSizeStruct.inGal.All.clMode=zeros(1,len);
cellSizeStruct.inGal.All.clQuant=zeros(4,len);
cellSizeStruct.inGal.All.clStd=zeros(1,len);
cellSizeStruct.inGal.All.clCount=zeros(1,len);

% cellSizeStruct.inGal.All.clTotalMean=0;
% cellSizeStruct.inGal.All.clTotalMean2=0;
% cellSizeStruct.inGal.All.clTotalSTD=0;
% cellSizeStruct.inGal.All.clTotalCount=0;


cellSizeStruct.inGal.SF=cellSizeStruct.inGal.All;
cellSizeStruct.inGal.NSF=cellSizeStruct.inGal.All;

cellSizeStruct.inCGM=cellSizeStruct.inGal;
cellSizeStruct.inSub=cellSizeStruct.inGal;
cellSizeStruct.inOut=cellSizeStruct.inGal;
cellSizeStruct.inAll=cellSizeStruct.inGal;

qus=[0.1 0.25 0.75 0.9];
cellSizeStruct.qants=qus;

%% loop over galaxies

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
        
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas',{'Coordinates','Density',...
        'Masses','StarFormationRate'});%'EnergyDissipation','Machnumber','Velocities',...
        
        
        if gas.count==0
            continue
        end
        
        
        cellSize=( double(gas.Masses) ./ double(gas.Density) ).^(1/3).*illUnits.lengthUnit;
        
        
        % find distance from galaxy center
        
        gas.newCoord = illustris.utils.centerObject(gas.Coordinates,subs.SubhaloPos(:,id+1));
        gasDist=sqrt( sum(double(gas.newCoord).^2,1));
        
        
        % find importat radii
        rmax=max(gasDist);
        rhalfStar=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('stars')+1,id+1)); % stellar half mass radius
        rhalfGas=double(subs.SubhaloHalfmassRadType(illustris.partTypeNum('gas')+1,id+1)); % gas half mass radius
        
        % identify Star-forming gas
        sfMask=gas.StarFormationRate>0;  %inpolygon(log10(nDensity),log10(gas.Temperature),polys.sf_polygonMask(:,1),polys.sf_polygonMask(:,2));
        
        
        for kk=1:5
            
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
                    
                case 5
                    % All gas in the Subfind
                    distMask=true(size(gasDist));
                    fld='inAll';
                    
            end
            
            
            for jj=1:3
                switch jj
                    case 1
                        % all cells
                        mask=distMask;
                        tag='All';
                    case 2
                        % all cells
                        mask=distMask & sfMask;
                        tag='SF';
                    case 3
                        % all cells
                        mask=distMask & ~sfMask;
                        tag='NSF';
                        
                end
                
                if any(mask)
                    cl=cellSize(mask);
                    
                    cellSizeStruct.(fld).(tag).clMean(id+1)=mean(cl);
                    cellSizeStruct.(fld).(tag).clMean2(id+1)=mean(cl.^2);
                    cellSizeStruct.(fld).(tag).clMedian(id+1)=median(cl);
                    cellSizeStruct.(fld).(tag).clQuant(:,id+1)=quantile(cl,qus);
                    cellSizeStruct.(fld).(tag).clStd(id+1)=calc_standardDev(cl);
                    cellSizeStruct.(fld).(tag).clCount(id+1)=sum(mask);
%                     cellSizeStruct.(fld).(tag).clTotalMean=...
%                         cellSizeStruct.(fld).(tag).clTotalMean+sum(cl);
%                     cellSizeStruct.(fld).(tag).clTotalMean2=...
%                         cellSizeStruct.(fld).(tag).clTotalMean2+sum(cl.^2);
%                     cellSizeStruct.(fld).(tag).clTotalCount=...
%                         cellSizeStruct.(fld).(tag).clTotalCount+sum(mask);
                    
                end
                
                
            end
            
        end
        
    end
    
end

% %% calculate the total mean and standard deviation
% 
% 
% for kk=1:5
%     
%     switch kk
%         case 1
%             fld='inGal';
%         case 2
%             fld='inCGM';
%         case 3
%             fld='inOut';
%         case 4
%             fld='inSub';
%         case 5
%             fld='inAll';
%     end
%     
%     for jj=1:3
%         switch jj
%             case 1
%                 tag='All';
%             case 2
%                 tag='SF';
%             case 3
%                 tag='NSF';
%         end
%         
%         cellSizeStruct.(fld).(tag).clTotalMean=...
%             cellSizeStruct.(fld).(tag).clTotalMean/cellSizeStruct.(fld).(tag).clTotalCount;
%         cellSizeStruct.(fld).(tag).clTotalMean2=...
%             cellSizeStruct.(fld).(tag).clTotalMean2/cellSizeStruct.(fld).(tag).clTotalCount;
%         cellSizeStruct.(fld).(tag).clTotalSTD=sqrt(...
%             cellSizeStruct.(fld).(tag).clTotalMean2-...
%             cellSizeStruct.(fld).(tag).clTotalMean^2);
%         
%     end
% end
% 

%% save to mat file

if DRACOFLAG
    %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fname=sprintf('gasCellSize_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'cellSizeStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
   
     
end
fprintf(' *** DONE!  *** \n');

%% save to postprocessed catalog

%global myPOSTPROCESSING
% cellSizeStruct.inGal.mask=floor(galMask);
% cellSizeStruct.inCGM.mask=floor(galMask);
% cellSizeStruct.inSub.mask=floor(galMask);
% illustris.utils.write_catalog( cellSizeStruct.inGal,99,'name','gasProps_inGal','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( cellSizeStruct.inCGM,99,'name','gasProps_inCGM','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( cellSizeStruct.inSub,99,'name','gasProps_inSub','folderName','gasProperties','verbose');
%



