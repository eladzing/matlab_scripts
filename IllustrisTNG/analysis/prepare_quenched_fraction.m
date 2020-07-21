%% analyze quenched fraction in clusters
for ii=1:2
    
    switch(ii)
        case 1
            sim='100';
        case 2
            sim='300';
    end
    
    
    
    snap=99; %z=0
    
    bp=illustris.set_env(sim);
    
    units; % load general unit structure in cgs.
    global illUnits
    global DEFAULT_MATFILE_DIR
    global LBox
    
    
    %load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat']);
     subs=illustris.groupcat.loadSubhalos(bp,snap);
     fofs=illustris.groupcat.loadHalos(bp,snap);
%     
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    %% set mask
    %%
    massThresh=10^9; % threshold for *stellar* mass
    
    %hostBins=[12 12.699 13 13.699 14 14.699 15 15.699];
    hostBins=12:0.5:15.5;
    
    galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'sats');
    
    
    ssfrThresh=1e-11;
    %% get data
    %%
    mass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galMask).*illUnits.massUnit;
    ssfr=illustris.utils.calc_ssfr(subs);
    ssfr=ssfr(galMask);
    hostMass=fofs.Group_M_Crit200(subs.SubhaloGrNr(galMask)+1).*illUnits.massUnit;
    %% set position in cluster
    %%
    satPos=subs.SubhaloPos(:,galMask);
    hostPos=fofs.GroupPos(:,subs.SubhaloGrNr(galMask)+1);
    
    dp0=abs(satPos-hostPos);
    dp=dp0;
    msk=dp>0.5*LBox;
    dp(msk)=LBox-dp(msk);
    
    hostRv=fofs.Group_R_Mean200(subs.SubhaloGrNr(galMask)+1);
    
    satRad=sqrt(sum((dp).^2,1))./hostRv;
    satRadProj=satRad.*generate_projectionFac(length(satRad));
    
    
    %% find fraction by distance host mass and stellar mass
    
    dr=0.25;
    rbin=0.25:dr:2.25;
    hmBin=hostBins;
    smBin=9:0.2:12;
    
    satRadbin=ceil(satRad./dr);
    satRadbin(satRadbin>length(rbin))=length(rbin);
    
    satRadbinProj=ceil(satRadProj./dr);
    satRadbinProj(satRadbinProj>length(rbin))=length(rbin);
    
    
    qbinT=zeros(length(smBin)-1,length(hmBin)-1,length(rbin));
    nbinT=qbinT;
    qbinTP=qbinT;
    nbinTP=qbinT;
    
    
    for k=1:length(smBin)-1
        ms=log10(mass)>=smBin(k) & log10(mass)<smBin(k+1);
        
        for j=1:length(hmBin)-1
            
            mh=log10(hostMass)>=hmBin(j) & log10(hostMass)<hmBin(j+1);
            
%             if j==1
%                 mh=log10(hostMass)<hmBin(j);
%             else
%                 mh=log10(hostMass)<hmBin(j) & log10(hostMass)>hmBin(j-1);
%             end
%             
            for i=1:length(rbin)
                m=satRadbin==i & mh & ms;
                mp=satRadbinProj==i & mh & ms;
                                
                qm=ssfr(m)<ssfrThresh;
                qbinT(k,j,i)=sum(qm);
                nbinT(k,j,i)=sum(m);
                
                qmP=ssfr(mp)<ssfrThresh;
                qbinTP(k,j,i)=sum(qmP);
                nbinTP(k,j,i)=sum(mp);
            end
            
        end
    end
    
    switch(ii)
        case 1
            qFracTNG.qbinT100=qbinT;
            qFracTNG.nbinT100=nbinT;
            qFracTNG.qbinProjT100=qbinTP;
            qFracTNG.nbinProjT100=nbinTP;
        case 2
            qFracTNG.qbinT300=qbinT;
            qFracTNG.nbinT300=nbinT;
            qFracTNG.qbinProjT300=qbinTP;
            qFracTNG.nbinProjT300=nbinTP;
    end
    
end
qFracTNG.hostBins=hostBins;
qFracTNG.stellarMassBins=smBin;
qFracTNG.radialBins=rbin;
%% save

fname=sprintf('quenchedFractions_r200M_TNG_snp%s',num2str(snap));
save([DEFAULT_MATFILE_DIR '/' fname],'qFracTNG','-v7.3')
