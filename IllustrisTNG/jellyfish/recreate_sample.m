for snap=[


illustris.utils.set_illUnits(snap);

global illUnits
if readFlag
    
    fprintf('Reading in data \n');
    
    % read in data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
end

global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/jf_objectTable_TNG50.mat']);
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/gasProperties_snp99_TNG300.mat')
    
% define stellar mass
starMassGal= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
starMassTot= double(subs.SubhaloMassType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass in subhalo 

% define gasMass 
gasMassGal= double(subs.SubhaloMassInRadType(illustris.partTypeNum('gas')+1,:).*illUnits.massUnit); % gas mass within 2*rhalf
gasMassTot= double(subs.SubhaloMassType(illustris.partTypeNum('gas')+1,:).*illUnits.massUnit); % gas mass in subhalo

%define host mass and radius 
r200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));
m200c=double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit);


% calculate distances between sub and host 
global LBox
rr=findDistance(double(subs.SubhaloPos), double(fofs.GroupPos(:,subsInfo.hostFof+1)),LBox);
rr1=findDistance(double(subs.SubhaloPos), double(fofs.GroupPos(:,subsInfo.hostFof+1)),1e10);

radPosition=rr./r200c;
radPosition1=rr1./r200c;
clear rr rr1 

%% generate masks 
galMassThresh=10^8.4;
hostMassThresh=10^11.5;
radLimit=1;
gasFracThresh=0.01;

% legal subhalos 
subMask=subs.SubhaloFlag;
fprintf('subMask pool = %i , %f \n',sum(subMask),sum(subMask)./length(subMask))

dmMask=subsInfo.DM10perc;
fprintf('dmMask pool = %i , %f \n',sum(dmMask),sum(dmMask)./length(dmMask))

% satellites 
satMask=~subsInfo.isCentral;
fprintf('satMask pool = %i , %f \n',sum(satMask),sum(satMask)./length(satMask))

%stellar mass
massMask=starMassGal>=galMassThresh;
fprintf('starMask pool = %i , %f \n',sum(massMask),sum(massMask)./length(massMask))

%host Mass
hostMask=m200c>=hostMassThresh;
fprintf('hostMask pool = %i , %f \n',sum(hostMask),sum(hostMask)./length(hostMask))

% distance Mask
distMask=radPosition<=radLimit;
distMask1=radPosition1<=radLimit;
fprintf('distMask pool = %i , %f \n',sum(distMask),sum(distMask)./length(distMask))
fprintf('distMask1 pool = %i , %f \n',sum(distMask),sum(distMask)./length(distMask))

% gasFrac mask 
gasMaskGal=(gasMassGal./starMassGal)>=gasFracThresh;
gasMaskTot=(gasMassTot./starMassTot)>=gasFracThresh;
gasMaskTry=(gasMassTot./starMassGal)>=gasFracThresh;
fprintf('gasMaskGal pool = %i , %f \n',sum(gasMaskGal),sum(gasMaskGal)./length(gasMaskGal))
fprintf('gasMaskTot pool = %i , %f \n',sum(gasMaskTot),sum(gasMaskTot)./length(gasMaskTot))
fprintf('gasMaskTry pool = %i , %f \n',sum(gasMaskTry),sum(gasMaskTry)./length(gasMaskTry))

%% combined sample 
mask=  satMask & massMask & hostMask & distMask & gasMaskGal & dmMask;
fprintf('total Mask pool = %i , %f \n',sum(mask),sum(mask)./length(mask))

mask1=  satMask & massMask & hostMask & distMask1 & gasMaskGal & dmMask;
fprintf('total Mask pool (no periodic cond.)= %i , %f \n',sum(mask1),sum(mask1)./length(mask1))

mask2=  satMask & massMask & hostMask & distMask1 & gasMaskGal & dmMask & subMask;
fprintf('total Mask pool (no periodic cond.)= %i , %f \n',sum(mask2),sum(mask2)./length(mask2))

fprintf('\n')

%% Check to see object by object if they exist 

sid=objectTable.subfind(objectTable.snap==snap); % subfind id's of objects in the table (zoo sample)
sidF=find(mask)'-1; %  % subfind id's of objects in larger sample 
sidF1=find(mask1)'-1;

found=false(size(sid));
for i=1:length(sid)
    found(i)=any(sidF==sid(i));
end

fprintf('Of %i objects in zoo sample, %i are in the current TNG sample (w/bound). \n \n',length(found),sum(found))


found1=false(size(sid));
for i=1:length(sid)
    found1(i)=any(sidF1==sid(i));
end

fprintf('Of %i objects in zoo sample, %i are in the current TNG sample (w/o bound). \n \n ',length(found1),sum(found1))


%% compile a list of objects not in the zoo sample. 

sidExtra=[];
for i=1:length(sidF)
    if ~any(sidF(i)==sid)
        sidExtra(end+1)=sidF(i);
    end
end

fprintf('Of %i objects in current sample (w/bound), %i are NOT in zoo sample . \n',length(sidF),length(sidExtra))



sidExtra1=[];
for i=1:length(sidF1)
    if ~any(sidF1(i)==sid)
        sidExtra1(end+1)=sidF1(i);
    end
end

fprintf('Of %i objects in current sample (w/o bound), %i are NOT in zoo sample . \n',length(sidF1),length(sidExtra1))

