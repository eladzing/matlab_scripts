%% This script generates a satellite catalog 

% Based on selection of random vantage point in the simulation and
% constructing everything compared to it.

% assume that we identify centrals ahead of time  

% sdss fiber effect - if more than galaxy is found within the aperature
% (distance dependent, only the most massive is recorded and the other is
% discarded 


% velocity dispersion is definded to be Vvir*factor 

% for a given central we find the los and projected distance of all other
% objects 
% 
% lucky centrals - consider 2 centrals, with the more massive one farther away than
% the other. The less massive may be beyond the excision aperature of the
% massive one and will not be removed. But the more massive may be in the
% excision aperature of the less massive one. For now we count the less
% massive one which should be removed - if the number is high we can remove
% them 

%% perliminary - load all necessary stuff
snap=99;
illustris.utils.set_illUnits(snap);
global illUnits
global cosmoStruct
global LBox
global DEFAULT_MATFILE_DIR
global simDisplayName;
 
zred=illustris.utils.snap2redshift(snap);

% set no. of vantage points (randomly selected)
Npoints=50;

% decide whether or not to include the sdss fiber effect 
if ~exist('fiberFlag','var') 
    error('Must define fiber flag')
end


if readFlag
    units;
    
    fprintf('Reading in data \n');
    
    % read in catalog data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    
end


%% set masks and central / satellite samples
virType='mean';

if setupFlag
     % define stellar mass
    galMass= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf

    massThresh=10^9;
    hostThresh=10^11;
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
    
    % define masks
    hostMask=m200>=hostThresh;  % cut hosts by M_200,c
    
    cenMask=centralMask & hostMask; 
    mm=centralMask & ~cenMask; % add centrals from smaller hosts to possible satellite pool
    satMask(mm)=true;
    clear mm 
    sampleMask= cenMask | satMask;
    
    
    %% use only relevant galaxies - reduces size of catalogs 
    galIDs=find(sampleMask)-1; % total gaalxy pool 
    galMass=galMass(sampleMask);
    r200=r200(sampleMask); % the host virial value for each gal 
    m200=m200(sampleMask);
    
    cenMask=cenMask(sampleMask);
    satMask=satMask(sampleMask);
    
    galPos=double(subs.SubhaloPos(:,sampleMask)); %3D position of galaxies 
    hostFof=subsInfo.hostFof(sampleMask); % host FOF ID for each galaxy 
    
    vsat=double(subs.SubhaloVel(:,sampleMask)).*illustris.utils.velocityFactor(snap,'sub'); % 3d velocity for each galaxy 
    vhost=double(fofs.GroupVel(:,subsInfo.hostFof(sampleMask)+1)).*illustris.utils.velocityFactor(snap,'host');% 3d velocity of host for each galaxy 
    
    %% rank by mass
    [galMass,ix]=sort(galMass,'descend');
    r200=r200(ix);
    m200=m200(ix);
    
    cenMask=cenMask(ix);
    satMask=satMask(ix);
    
    galPos=galPos(:,ix);
    hostFof=hostFof(ix);
    vsat=vsat(:,ix);
    vhost=vhost(:,ix);
end

%% set fiber excision paramters
% translate 33'' to kpc
zz=linspace(0,0.1,100); % redshift array 
ang=33/3600*pi/180;
angPhys=ang.*angular_distance(zz,'cosmo',cosmoStruct)*1000; % size of aperture (in physical kpc) as a function of redshift

% use comoving distance to translate distance to redshift
zlen=comoving_distance(zz);  
% interplation on the  angPhys vs. zlen relation will give the fiber
% aperature for a given object separation 

% alternate method 
% ll=linspace(0,300,1000); % in Mpc
% zlen2=ll.*(cosmoStruct.hub*100).*Units.km./Units.cspeed;


%% set Yang model stuff 
sigmaFac=1;%/sqrt(2);
radFac=2; % 
bValue=1;  % threshold for Yang model 10 corresponds to ~rvir, 1 corresponds to ~2 Rvir 

v200=sqrt(Units.GG.*m200./r200); % in km/sec;
vSigma=v200.*sigmaFac; % how do you define the velocity dispersion 

cv=cvir_Mvir(m200,zred,'hub',cosmoStruct.hub); % concentration to mvir relation 
delta=200/3.*cv.^3./(log(1+cv)-cv./(1+cv)); % NFW definition of overdensity 

rs=r200./cv; %in kpc - NFW scale radius 
bigSigmaFac=2.*(rs.*1e-3).*delta; 

%% select vantage point(s)
center=rand(3,Npoints).*LBox; %in kpc/h 

%% loop over vantage points

for k=1:Npoints
    
    cMask=cenMask;
    sMask=satMask;
    
    luckyCentral=[]; % centrals which should be excised by fiber but arent
    
    %new coordinates in kpc - center everything around vantage point
    newCoords=illustris.utils.centerObject(galPos,center(:,k)).*illUnits.lengthUnit;
    
    
    %% run over centrals 
    stp=10;
    prc=10;
    cntr=0;
    
    scnt=1;
    for j=1:length(cMask)
        
        % skip gals not on the central list 
        if ~cMask(j) 
            continue
        end
        
        % print counter
        cntr=cntr+1;
        if floor(cntr./sum(cMask).*100)>=prc
            fprintf('completed %s %% of centrals \n',num2str(prc));
            prc=prc+stp;
        end
        
        
        
        %% for each central - set projected and los distances of other galaxies
        % center everone around central
        rCen=newCoords(:,j); % position vector of central w.r.t vantage point
        distCen=sqrt(sum(rCen.^2));
        rCenHat=rCen./distCen;
        
        % position vector of all other galaxies w.r.t the central
        satCoords=illustris.utils.centerObject(galPos,galPos(:,j))...
            .*illUnits.lengthUnit;
        satDist=sqrt(sum(satCoords.^2,1));
        
        % calculate projected and los distnce
        losComp=satCoords(1,:).*rCenHat(1) + satCoords(2,:).*rCenHat(2) + satCoords(3,:).*rCenHat(3);
        projComp=sqrt(satDist.^2-losComp.^2);
        
        
        %% for each central, excise all galaxies which are too close - fiber effect
        if fiberFlag
            fiberAper=interp1(zlen,angPhys,distCen./1000);
            fiberMask= projComp<= fiberAper & projComp>0;
        else
            fiberMask=false(size(projComp));
        end
        
        sMask(fiberMask)=false;  % remove these galaxies from sample
        cMask(fiberMask)=false;  % remove these galaxies from sample
        
        % count centrals whish should be removed but aren't
        if any(galMass(fiberMask)>galMass(j))
            luckyCentral(end+1)=j;
            % if we wish to reomve them: 
            % cMask(j)=false;
            % continue
        end
        
        %% generate Yang score for all satellites
        bigSigma=bigSigmaFac.*auxSigma(projComp./rs(j));
        distMask=projComp<=radFac.*r200(j); % projected distance cut
        
        % Hubble flow
        vhub=(cosmoStruct.hub.*100).*(losComp./1000);
        % set to zero if the satellite is part of the
        % group:        
        % actual members of the host fof
        memberMask=hostFof==hostFof(j);
        vhub(memberMask)=0;
        
        %los velocity
        vpec=vsat-vhost(j);
        vlos=vpec(1,:).*rCenHat(1) + vpec(2,:).*rCenHat(2) + vpec(3,:).*rCenHat(3);
        
        vv=vhub+vlos;
        pm=pmax(vv,vSigma(j),zred);
        pscore=(cosmoStruct.hub.*100).*bigSigma.*pm;
        
        % keep only satellites above threshold
        satMaskYang=sMask & pscore>bValue & ~fiberMask;
        
        if sum(satMaskYang)>0
            inds=scnt:scnt+sum(satMaskYang)-1;
            
            scnt=inds(end)+1;
            
            satStructY(k).satID(inds)=galIDs(satMaskYang);% find(satMaskYang)-1;
            satStructY(k).hostID(inds)=hostFof(j);
            satStructY(k).cenID(inds)=galIDs(j);
            satStructY(k).distProj(inds)=projComp(satMaskYang);
            satStructY(k).pscore(inds)=pscore(satMaskYang);
            
            
            
        end
        
    end
    
    satStructY(k).luckyCentralInd=luckyCentral;
    satStructY(k).luckyCentralMass=galMass(luckyCentral);
    
    %% reomve double counts based on maximal score.
    
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
    
            
            satStructY(k).satID=satStructY(k).satID(listMask);
            satStructY(k).hostID=satStructY(k).hostID(listMask);
            satStructY(k).cenID=satStructY(k).cenID(listMask);
            satStructY(k).distProj=satStructY(k).distProj(listMask);
            satStructY(k).pscore=satStructY(k).pscore(listMask);
    
fprintf('completed %i vantage points out of %i \n',k,Npoints);
end

if fiberFlag
    catName=[DEFAULT_MATFILE_DIR '/yangSatSample_Vpoint_sig1_' virType '_fiber_' simDisplayName];
else
    catName=[DEFAULT_MATFILE_DIR '/yangSatSample_Vpoint_sig1_' virType '_' simDisplayName];
end

save(catName,'satStructY');

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

