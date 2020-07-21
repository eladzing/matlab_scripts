%% this is an updated version of generate_coolingTimes_gals
% run over all subs and extract radius of farthest gas cell.

%% set framework

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
    
end


massThresh=10^9; % threshold for *stellar* mass

%% load subsinfo
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.
%% select galaxies

massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf

% this mask selects all galaxies with dm component & stars, above *stellar* mass limit whose host has virials parameters
galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
radiiStruct.galMass=massAllGals;
radiiStruct.galMask=galMask;

r200c=fofs.Group_R_Crit200(subsInfo.hostFof+1);

%% generate values
fprintf(' *** Running over %s Galaxies *** \n',num2str(sum(galMask)));

step=5;
stepNext=5;
len=double(subs.count);

%% initialize to zero

radiiStruct.rStar=zeros(1,len);
radiiStruct.rhalfStar=zeros(1,len);
radiiStruct.rhalfGas=zeros(1,len);
radiiStruct.rMaxGas=zeros(1,len);
radiiStruct.r90Gas=zeros(1,len);
radiiStruct.r95Gas=zeros(1,len);
%radiiStruct.massFunction(2,len);


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
        
        gas=illustris.snapshot.loadSubhalo(bp, snap, id, 'gas',{'Coordinates','Masses'});
        
        if gas.count==0
            continue
        end
        
        
        mass=double(gas.Masses.*illUnits.massUnit); %in Solarmass
        
        totalMass=sum(mass);
        
        
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
        
        if gas.count>2
        
        [sDist, sInd]=sort(gasDist);
        csMass=cumsum(mass(sInd))./totalMass;
        
        r90=interp1(csMass,sDist,0.9);
        r95=interp1(csMass,sDist,0.95);
        else
            r90=0;
            r95=0;
        end
        
        radiiStruct.rMaxGas(id+1)=rmax.*illUnits.lengthUnit;
        radiiStruct.r90Gas(id+1)=r90.*illUnits.lengthUnit;
        radiiStruct.r95Gas(id+1)=r95.*illUnits.lengthUnit;
        radiiStruct.rhalfStar(id+1)=rhalfStar.*illUnits.lengthUnit;% in kpc
        radiiStruct.rStar(id+1)=2.0*rhalfStar.*illUnits.lengthUnit;% in kpc
        radiiStruct.rhalfGas(id+1)=rhalfGas.*illUnits.lengthUnit; % in kpc
        radiiStruct.gasCount(id+1)=gas.count;
       
    end
end


%% save to mat file
global DRACOFLAG
if DRACOFLAG
    global DEFAULT_MATFILE_DIR
    global simDisplayName
    %fname=sprintf('gasProperties_z%s_%s',num2str(illustris.utils.get_zred(snap)),simDisplayName);
    fname=sprintf('gasRadii_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'radiiStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    
end
fprintf(' *** DONE!  *** \n');

%% save to postprocessed catalog

%global myPOSTPROCESSING
% radiiStruct.inGal.mask=floor(galMask);
% radiiStruct.inCGM.mask=floor(galMask);
% radiiStruct.inSub.mask=floor(galMask);
% illustris.utils.write_catalog( radiiStruct.inGal,99,'name','gasProps_inGal','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( radiiStruct.inCGM,99,'name','gasProps_inCGM','folderName','gasProperties','verbose');
% illustris.utils.write_catalog( radiiStruct.inSub,99,'name','gasProps_inSub','folderName','gasProperties','verbose');
%



