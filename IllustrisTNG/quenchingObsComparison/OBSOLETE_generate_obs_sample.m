snap=99;
illustris.utils.set_illUnits(snap);
global illUnits
global cosmoStruct
units;
%%  read stuff

massThresh=10^9;
hostThresh=10^11;

if readFlag
    
    fprintf('Reading in data \n');
    
    % read in data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
end
if setupFlag
    % define stellar mass
    massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
    
    sampleMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
    centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
    
    %define host mass and radius
    r200=double(fofs.Group_R_Mean200(subsInfo.hostFof+1)).*illUnits.lengthUnit;
    m200=double(fofs.Group_M_Mean200(subsInfo.hostFof+1).*illUnits.massUnit);
    
    r200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));
    m200c=double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit);
    
    v200=sqrt(Units.GG.*m200./r200);
    v200c=sqrt(Units.GG.*m200c./(r200c.*illUnits.lengthUnit));
    sigma=v200c;%./sqrt(2);
    
    aexp=illustris.utils.snap2aexpn(snap);
end


%% run over all centrals

radFac=2;

carts=1:3;
cnt(1:3)=1;

hostMask=m200c>=hostThresh;

cenMask=centralMask & hostMask;
sampleMask=sampleMask & ~cenMask;

stp=10;
prc=10;
cntr=0;

for j=1:length(centralMask)
    
    % skip
    if ~cenMask(j)
        continue
    end
    cntr=cntr+1;
    
    if floor(cntr./sum(cenMask).*100)>=prc
        fprintf('completed %s %% of centrals \n',num2str(prc));
        prc=prc+stp;
    end
    
    hostInd=subsInfo.hostFof(j)+1;
    
    % center everyone around the central (takes care of boundary
    % conditions)
    fofCen=double(fofs.GroupPos(:,hostInd));
    newCoords=illustris.utils.centerObject(double(subs.SubhaloPos),fofCen);
    
    memberMask=subsInfo.hostFof==subsInfo.hostFof(j);

    
    for k=1:3 % generate for each cartesian direction
        
        kk=carts(carts~=k);
        
        % transverse distance (normalized)
        rr=sqrt(newCoords(kk(1),:).^2+newCoords(kk(2),:).^2).*illUnits.lengthUnit;
        
        %l.o.s length
        zlen=double(newCoords(k,:).*illUnits.lengthUnit);
        
        % Hubble flow
        vhub=(cosmoStruct.hub.*100).*(zlen./1000);
        vhub(memberMask)=0;
        
        %los velocity
        vlos=(double(subs.SubhaloVel(k,:))-double(fofs.GroupVel(k,hostInd))./aexp);
        
        velMask=abs((vhub+vlos))<=sigma(j);
        distMask= rr<=(radFac.*r200(j));
        
        
        satMask=sampleMask & distMask & velMask;
        satMask(j)=false;
        
        if sum(satMask)>0
            inds=cnt(k):cnt(k)+sum(satMask)-1;
            
            cnt(k)=inds(end)+1;
            
            satStruct(k).satID(inds)=find(satMask)-1;
            satStruct(k).hostID(inds)=hostInd-1;
            satStruct(k).cenID(inds)=j-1;
            satStruct(k).distProj(inds)=rr(satMask);
            
        end
    end
end

%% remove double counting 
for k=1:3
    listMask=false(size(satStruct(k).satID));
    
    for j=1:length(satStruct(k).satID)
        
        indlist=find(satStruct(k).satID==satStruct(k).satID(j));
        if length(indlist)==1
            listMask(indlist)=true;
        else
            ps=satStruct(k).pscore(indlist);
            [~,ix]=max(ps);
            listMask(indlist(ix))=true;
            
        end
        
    end
    
    satStruct(k).listMask=listMask;
end



%% save information as table and in mat file
for k=1:3
    msk=satStruct(k).listMask;
    satTab=table(int64(satStruct(k).satID(msk)'),satStruct(k).hostID(msk)',satStruct(k).cenID(msk)',satStruct(k).distProj(msk)',...
        'variableNames',{'satID','hostID','centralID','Projected_Distance'});
    switch k
        case 1
            satTabX=satTab;
        case 2
            satTabY=satTab;
        case 3
            satTabZ=satTab;
    end
end


global DEFAULT_MATFILE_DIR
global simDisplayName;

catName=[DEFAULT_MATFILE_DIR '/SatSample_vCut_sig1_' simDisplayName];
save(catName,'satTabX','satTabY','satTabZ')
fprintf('sat catalog written to %s \n',catName);

