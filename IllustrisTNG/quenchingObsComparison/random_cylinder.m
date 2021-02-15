snap=99;
illustris.utils.set_illUnits(snap);

%%  read stuff 

if readFlag
    
    fprintf('Reading in data \n');
    
    % read in data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/gasProperties_snp99_TNG300.mat')
    
    % define stellar mass
    massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
    
    sampleMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);
    centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
    
    %define host mass and radius 
    r200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));
    m200c=log10(double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit));

    ssfr=double(illustris.utils.calc_ssfr(subs,'base',0));

end
    


%% run over all centrals above host mass cut 
binEdges=0.1:0.2:2;

radFac=2;
hostRange=11:15;
hostRange=10.^cat(2,hostRange',hostRange'+1);

% zThresh=  This needs to be defined base don Obs cuts.   

carts=1:3;
for i=1:length(hostRange)
    hostMask=m200c>=hostRange(i,1) & m200c<=hostRange(i,2);
    
    loopMask=centralMask & hostMask;
    
    for j=1:length(centralMask)
        
        % skip 
        if loopMask(j)
            continue
        end
        
        newCoords=illustris.utils.centerObject(subs.SubhaloPos,subs.SubhaloPos(:,j));
        
        
        
        for k=1:3 % generate for each cartesian direction 
            
            kk=carts(carts~=k);r
            rr=sqrt(newCoords(kk(1),:).^2+newCoords(kk(2),:).^2)./r200c;
            zlen=newCoords(k,:);
                    
            satMask=sampleMask & rr<=radFac & abs(zLen)<zThresh(j); 
            

%% center all subs around central in question 

%% generate random 
% phiRange=[0 2*pi];
% psiRange=phiRange;
% zRange=[-1 1];
% %thetRange=[-pi/2 pi/2];
% 
% res.phi=phiRange(1)+diff(phiRange).*rand(nn,1);
% res.psi=psiRange(1)+diff(psiRange).*rand(nn,1);
% zz=zRange(1)+diff(zRange).*rand(nn,1);
% 
% res.theta=acos(zz);

%% 

