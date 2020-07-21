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
    %global LBox
    
    
    %load([DEFAULT_MATFILE_DIR '/fofs_subs_TNG100_z0.mat']);
    subs=illustris.groupcat.loadSubhalos(bp,snap);
    fofs=illustris.groupcat.loadHalos(bp,snap);
%     
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    %% set mask
    %%
    massThresh=10^9; % threshold for *stellar* mass
    
    %hostBins=[12 12.699 13 13.699 14 14.699 15 15.699];
    hostBins=[10.9 11.5:0.5:15.5];
    
    spbMask = mk_splashback_mask('time',5,'both',0.1);
    
    gMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
    
    galMask=gMask & ~spbMask;
       
    
    ssfrThresh=1e-11;
    %% get data
    %%
    mass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,galMask).*illUnits.massUnit;
    ssfr=illustris.utils.calc_ssfr(subs);
    ssfr=ssfr(galMask);
    hostMass=fofs.Group_M_Crit200(subs.SubhaloGrNr(galMask)+1).*illUnits.massUnit;
    %% set position in cluster
    %%

    hostRv=fofs.Group_R_Crit200(subs.SubhaloGrNr(galMask)+1);
    
  
    
    %% find fraction by distance host mass and stellar mass
    

    hmBin=hostBins;
    smBin=9:0.2:12;
    

    
    qbinT=zeros(length(smBin)-1,length(hmBin)-1);
    nbinT=qbinT;
    
    
    
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
             m= mh & ms;
                
                                
                qm=ssfr(m)<ssfrThresh;
                qbinT(k,j)=sum(qm);
                nbinT(k,j)=sum(m);


   
        end
    end
    
    switch(ii)
        case 1
            qFracCenTNG.qbinT100=qbinT;
            qFracCenTNG.nbinT100=nbinT;
            
        case 2
            qFracCenTNG.qbinT300=qbinT;
            qFracCenTNG.nbinT300=nbinT;
            
    end
    
end
qFracCenTNG.hostBins=hostBins;
qFracCenTNG.stellarMassBins=smBin;
%% save

fname=sprintf('quenchedFractionsCentrals_TNG_z%s',num2str(illustris.utils.get_zred(snap)));
save([DEFAULT_MATFILE_DIR '/' fname],'qFracCenTNG','-v7.3')
