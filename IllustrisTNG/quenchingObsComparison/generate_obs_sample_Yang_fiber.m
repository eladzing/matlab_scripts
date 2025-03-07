b%% generates mock catalogs of projected satellite-halo catalogs  using three different methods: 
% 1) recreating the Yang et al. 2005 method
% 2) projected radius at 2 X Rvir & velocity cut within 0.5 vvir
% 3) projected radius at 2 X Rvir & within redshift cylinder set by vvir 

%% perliminaries 

snap=99;
illustris.utils.set_illUnits(snap);
global illUnits
global cosmoStruct
global LBox
units;


sigmaFac=1;%/sqrt(2);
radFac=2;


virType='mean';
%%  read stuff

massThresh=10^9;
hostThresh=10^11;

bValue=1;  % threshold for Yang model .

zred=illustris.utils.snap2redshift(99);

if readFlag
    
    fprintf('Reading in data \n');
    
    % read in data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
end

if setupFlag
    % define stellar mass
    massAllGals= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf
    
    %sampleMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
    centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
    satMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'sats');
    
    %define host mass and radius
    r200m=double(fofs.Group_R_Mean200(subsInfo.hostFof+1)).*illUnits.lengthUnit;
    m200m=double(fofs.Group_M_Mean200(subsInfo.hostFof+1).*illUnits.massUnit);

    r200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1)).*illUnits.lengthUnit;
    m200c=double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit);
        
    aexp=illustris.utils.snap2aexpn(snap);

    switch virType
        case 'mean'
            r200=r200m;
            m200=m200m;
            
        case 'crit'
            r200=r200c;
            m200=m200c;
    end
    
    
    
    hostMask=m200c>=hostThresh;  % cut hosts by M_200,c 

    cenMask=centralMask & hostMask;
    mm=centralMask & ~cenMask;
    satMask(mm)=true;
    sampleMask= cenMask | satMask;
    
    
end 

if fiberFlag
   % translate 33'' to kpc
   zz=linspace(0,0.075,100);
   ang=33/3600*pi/180;
   angPhys=ang.*angular_distance(zz,'cosmo',cosmoStruct)*1000; % in kpc

   % use hubble law to translate distance to redshift
   ll=linspace(0,300,1000); % in Mpc
   zl=ll.*(cosmoStruct.hub*100).*Units.km./Units.cspeed;
   
   kk=[2 3 ; 1 3 ; 1 2];
   for k=1:3  % loop over projections
       
       % projected distance
       pos(1,:)=double(subs.SubhaloPos(kk(k,1),:)./LBox);
       pos(2,:)=double(subs.SubhaloPos(kk(k,2),:)./LBox);
       pos3=double(subs.SubhaloPos(k,:).*illUnits.lengthUnit)./1000;
       %zlen=interp1(ll,zl,pos3,'linear');
       
       % constant aperature for z=0.03
       fibarPar=0.03;
       lp=interp1(zz,angPhys,fibarPar,'linear')./(LBox.*illUnits.lengthUnit.*(1/(1+fibarPar)));
       
       fiber(k)=sdssFiberSelect(pos,massAllGals,lp,sampleMask);
       
   end
else
    fiber(1).outMask=true(size(massAllGals));
    fiber(2).outMask=true(size(massAllGals));
    fiber(3).outMask=true(size(massAllGals));
end




v200=sqrt(Units.GG.*m200./r200); % in km/sec;
vSigma=v200.*sigmaFac;


% yang stuff
cv=cvir_Mvir(m200,zred,'hub',cosmoStruct.hub);
delta=200/3.*cv.^3./(log(1+cv)-cv./(1+cv));

rs=r200./cv; %in kpc
bigSigmaFac=2.*(rs.*1e-3).*delta;

% redshift cylinder
deltaZ=vSigma./(Units.cspeed./Units.km);




%% run over all centrals

carts=1:3;
cntY(1:3)=1;
cntZ(1:3)=1;
cntV(1:3)=1;



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
    
    
    % host velocity dispersion
    %sigma(1)=v200(j);
    %sigma(2)=v200(hostInd)./sqrt(2);
    %sigma(3)=397.9.*(m200(hostInd)./(1e14./cosmoStruct.hub)).^0.3124;
    
    % actual members of the host fof
    memberMask=subsInfo.hostFof==subsInfo.hostFof(j);
    
    for k=1:3 % generate for each cartesian direction
        
        if ~fiber(k).outMask(j) % skip if central is removed by fiber 
            continue
        end
        
        kk=carts(carts~=k);
        
        %% transverse distance
        rr=sqrt(newCoords(kk(1),:).^2+newCoords(kk(2),:).^2).*illUnits.lengthUnit;
        
        % yang 
        bigSigma=bigSigmaFac.*auxSigma(rr./rs(j));
        
        distMask=rr<=radFac.*r200(j); % projected distance cut 
        
        %% velocity
        %l.o.s length
        zlen=double(newCoords(k,:).*illUnits.lengthUnit);
        
        % redshift cylinder 
        losLength=comoving_distance(deltaZ(j),'cosmo',cosmoStruct).*1000; % in kpc
        losLengthMask=abs(zlen)<=0.5.*losLength;
        
        % Hubble flow
        vhub=(cosmoStruct.hub.*100).*(zlen./1000);
        % perhaps we should set to zero if the satellite is part of the
        % group.
        vhub(memberMask)=0;
               
        %los velocity
        vsat=double(subs.SubhaloVel(k,:)).*illustris.utils.velocityFactor(snap,'sub')
        vlos=(-double(fofs.GroupVel(k,hostInd))./aexp);
        vv=vhub+vlos;
        velMask=abs(vv)<=vSigma(j);
        
        % yang pscore 
         pm=pmax(vv,vSigma(j),zred);
         pscore=(cosmoStruct.hub.*100).*bigSigma.*pm;
        
         % Masks 
        satMaskYang=satMask & pscore>bValue & fiber(k).outMask;
        %satMaskYang(j)=false;  % don'tcount central as satellite. 
        
        satMaskVel=satMask & distMask & velMask & fiber(k).outMask;
        %satMaskVel(j)=false;
        
        satMaskZred=satMask & distMask & losLengthMask & fiber(k).outMask;
        %satMaskZred(j)=false;
        
        if sum(satMaskYang)>0
            inds=cntY(k):cntY(k)+sum(satMaskYang)-1;
            
            cntY(k)=inds(end)+1;
            
            satStructY(k).satID(inds)=find(satMaskYang)-1;
            satStructY(k).hostID(inds)=hostInd-1;
            satStructY(k).cenID(inds)=j-1;
            satStructY(k).distProj(inds)=rr(satMaskYang);
            satStructY(k).pscore(inds)=pscore(satMaskYang);
        end
        
         if sum(satMaskVel)>0
            inds=cntV(k):cntV(k)+sum(satMaskVel)-1;
            
            cntV(k)=inds(end)+1;
            
            satStructV(k).satID(inds)=find(satMaskVel)-1;
            satStructV(k).hostID(inds)=hostInd-1;
            satStructV(k).cenID(inds)=j-1;
            satStructV(k).distProj(inds)=rr(satMaskVel);
            satStructV(k).vv=abs(vv(inds));
         end
        
         if sum(satMaskZred)>0
            inds=cntZ(k):cntZ(k)+sum(satMaskZred)-1;
            
            cntZ(k)=inds(end)+1;
            
            satStructZ(k).satID(inds)=find(satMaskZred)-1;
            satStructZ(k).hostID(inds)=hostInd-1;
            satStructZ(k).cenID(inds)=j-1;
            satStructZ(k).distProj(inds)=rr(satMaskZred);
            satStructZ(k).zlen=abs(zlen(inds));
        end
        
    end
end

%% remove double counting Yang
for k=1:3
    listMask=false(size(satStructY(k).satID));
    skipMask=listMask;
    
    for j=1:length(satStructY(k).satID)
        
        if skipMask(j)
            continue
        end
        
        indlist=find(satStructY(k).satID==satStructY(k).satID(j));
        skipMask(indlist)=true;
        if length(indlist)==1
            listMask(indlist)=true;
        else
            ps=satStructY(k).pscore(indlist);
            [~,ix]=max(ps);
            listMask(indlist(ix))=true;
            
        end
        
    end
    
    satStructY(k).listMask=listMask;
end


%% remove double counting zred
for k=1:3
    listMask=false(size(satStructZ(k).satID));
    skipMask=listMask;
    
    for j=1:length(satStructZ(k).satID)
        
        if skipMask(j)
            continue
        end
        
        indlist=find(satStructZ(k).satID==satStructZ(k).satID(j));
        skipMask(indlist)=true;
        if length(indlist)==1
            listMask(indlist)=true;
        else
            rp=sqrt(satStructZ(k).distProj(indlist).^2+satStructZ(k).zlen.^2);
            [~,ix]=min(rp);
            listMask(indlist(ix))=true;
            
        end
        
    end
    
    satStructZ(k).listMask=listMask;
end

%% remove double counting vel
for k=1:3
    listMask=false(size(satStructV(k).satID));
    skipMask=listMask;
    
    for j=1:length(satStructV(k).satID)
        
        if skipMask(j)
            continue
        end
        
        indlist=find(satStructV(k).satID==satStructV(k).satID(j));
        skipMask(indlist)=true;
        if length(indlist)==1
            listMask(indlist)=true;
        else
            rp=sqrt(satStructV(k).distProj(indlist).^2+satStructV(k).vv.^2);
            [~,ix]=min(rp);
            listMask(indlist(ix))=true;
            
        end
        
    end
    
    satStructV(k).listMask=listMask;
end


%% save information as table and in mat file Yang
global DEFAULT_MATFILE_DIR
global simDisplayName;


for k=1:3
    msk=satStructY(k).listMask;
    satTab=table(int64(satStructY(k).satID(msk)'),satStructY(k).hostID(msk)',satStructY(k).cenID(msk)',satStructY(k).distProj(msk)',satStructY(k).pscore(msk)',...
        'variableNames',{'satID','hostID','centralID','Projected_Distance','pscore'});
    switch k
        case 1
            satTabX=satTab;
        case 2
            satTabY=satTab;
        case 3
            satTabZ=satTab;
    end
end

catName=[DEFAULT_MATFILE_DIR '/yangSatSample_fiber_sig1_' virType '_' simDisplayName];
save(catName,'satTabX','satTabY','satTabZ');

fprintf('sat catalog written to %s \n',catName);


%% save information as table and in mat file Vel

for k=1:3
    msk=satStructV(k).listMask;
    satTab=table(int64(satStructV(k).satID(msk)'),satStructV(k).hostID(msk)',satStructV(k).cenID(msk)',satStructV(k).distProj(msk)',...
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

catName=[DEFAULT_MATFILE_DIR '/velSatSample_fiber_sig1_' virType '_' simDisplayName];
save(catName,'satTabX','satTabY','satTabZ');

fprintf('sat catalog written to %s \n',catName);

%% save information as table and in mat file Vel

for k=1:3 
    msk=satStructZ(k).listMask;
    satTab=table(int64(satStructZ(k).satID(msk)'),satStructZ(k).hostID(msk)',satStructZ(k).cenID(msk)',satStructZ(k).distProj(msk)',...
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

catName=[DEFAULT_MATFILE_DIR '/zredSatSample_fiber_sig1_' virType '_' simDisplayName];
save(catName,'satTabX','satTabY','satTabZ');

fprintf('sat catalog written to %s \n',catName);





%% auxilary functions

function res=auxSigma(xx)

res=zeros(size(xx));

msk=xx>1;

xx2=xx(msk).^2-1;
res(msk)=(1 - atan(sqrt(xx2))./sqrt(xx2) )./(xx2);

msk=xx<1;
xx2=1-xx(msk).^2;
res(msk)=(log( (1+sqrt(xx2))./xx(msk) )./sqrt(xx2)-1)./(xx2);

msk=xx==1;
res(msk)=1/3;

end


function res=pmax(dv,sig,zgr)

%cc=Units.cspeed./Units.km;
cc=1;
fac=cc/sqrt(2*pi);
sigg=sig.*(1+zgr);
quot=-0.5.*(dv./sigg).^2;
res=fac.*exp(quot)./sigg;

end
