function ssfrMassProfs= generate_profile_from_vantagePoint_Catalog(satCat,fofs,subs,varargin)



global illUnits
% global DEFAULT_MATFILE_DIR
% global simDisplayName

massThresh=10^9;
virType='crit200';
snap=99;
%% parse arguments 
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'mass','massthresh','thresh','threshold'}
           i=i+1;
            massThresh=varargin{i};
        case {'crit','crit200','200crit'}
            virType='crit200';
        case {'mean','mean200','200mean'}
            virType='mean200';
        case {'500','crit500','500crit'}
            virType='crit500';
        case 'snap'
            i=i+1;
            snap=varargin{i};
            
            
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end

%% Perliminaries  
illustris.utils.set_illUnits(snap);
% read in data
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
%load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/gasProperties_snp99_TNG300.mat')

% define stellar mass
massAllGals= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf

%% calculate sSFR
ssfrBase=double(illustris.utils.calc_ssfr(subs,'base',0));

%% set mask
%satMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'sats');
%     centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
%


%% define virial parameters
switch lower(virType)
    case 'crit200'
        rvir=double(fofs.Group_R_Crit200(subsInfo.hostFof+1)).*illUnits.lengthUnit; % host rvir for each galaxy
        mvir=double(fofs.Group_M_Crit200.*illUnits.massUnit);                       % host mvir for each Fof
    case 'mean200'
        rvir=double(fofs.Group_R_Mean200(subsInfo.hostFof+1)).*illUnits.lengthUnit;
        mvir=double(fofs.Group_M_Mean200.*illUnits.massUnit); 
    case 'crit500'
        rvir=double(fofs.Group_R_Crit500(subsInfo.hostFof+1)).*illUnits.lengthUnit;
        mvir=double(fofs.Group_M_Crit500.*illUnits.massUnit); 
end
%     m200c=log10(double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit));


%load([DEFAULT_MATFILE_DIR '/yangSatSample_Vpoint_sig1_mean_' simDisplayName '.mat'])
%fname=[DEFAULT_MATFILE_DIR '/ssfr0_stellarMass_radProfiles_projected_yangSamplePoint_' simDisplayName];

illustris.utils.set_illUnits(snap);

%% bin the distances
binEdges=0.1:0.2:2;
%binInd=discretize(radPosition,binEdges);
%rbin=binEdges(1:end-1)+0.5.*diff(binEdges);
%% run with all galaxies - includes ssfr=0

% %bin by host mass;
% mask11=m200c>=11 & m200c<12;
% mask12=m200c>=12 & m200c<13;
% mask13=m200c>=13 & m200c<14;
% mask14=m200c>=14 & m200c<15;
% 
% % central masks
% cmask11=centralMask & m200c>=11 & m200c<12;
% cmask12=centralMask & m200c>=12 & m200c<13;
% cmask13=centralMask & m200c>=13 & m200c<14;
% cmask14=centralMask & m200c>=14 & m200c<15;

%% initialize
xMedian=zeros(length(satCat),length(binEdges)-1,4);
xMean=zeros(length(satCat),length(binEdges)-1,4);

smassMedian=zeros(length(satCat),length(binEdges)-1,4);
smassMean=zeros(length(satCat),length(binEdges)-1,4);

ssfrMedian=zeros(length(satCat),length(binEdges)-1,4);
ssfrMean=zeros(length(satCat),length(binEdges)-1,4);

for k=1:length(satCat)
    % enforce mass threshold
    satInd=satCat(k).satID+1;
    mass=massAllGals(satInd);
    massMask=mass>=massThresh;
    tempCat=mask_structure(satCat(k),massMask);
    
    
    satInd=tempCat.satID+1;
    mass=massAllGals(satInd);
    
    hostInd=tempCat.hostID+1;  
        
    ssfr=ssfrBase(satInd);
    radPosition=tempCat.distProj./rvir(satInd);
    hostMass=log10(mvir(hostInd));
    
    hmask=false(4,length(hostMass));        
    hmask(1,:)=hostMass>=11 & hostMass<12;
    hmask(2,:)=hostMass>=12 & hostMass<13;
    hmask(3,:)=hostMass>=13 & hostMass<14;
    hmask(4,:)=hostMass>=14 & hostMass<15;
    
    for j=1:4
        mask=squeeze(hmask(j,:));
        ssfrProf=mk_meanMedian_bin(radPosition(mask),(ssfr(mask)),'bins',binEdges);
        starMassProf=mk_meanMedian_bin(radPosition(mask),(mass(mask)),'bins',binEdges);
        
        ssfrAvg(j)=ssfrProf;
        starMass(j)=starMassProf;
    
    
        %% populate arrays for calculating global profile
      
        xMedian(k,:,j)=starMassProf.xMedian;
        xMean(k,:,j)=starMassProf.xMean;
        
        smassMedian(k,:,j)=starMassProf.yMedian;
        smassMean(k,:,j)=starMassProf.yMean;
       
       
        ssfrMedian(k,:,j)=ssfrProf.yMedian;
        ssfrMean(k,:,j)=ssfrProf.yMean;
    end
    
    
    
%     % centrals
%     ssfrAvgC(1)=mk_meanMedian_bin(radPosition(cmask11),(ssfr(cmask11)),'bins',[0 1]);
%     ssfrAvgC(2)=mk_meanMedian_bin(radPosition(cmask12),(ssfr(cmask12)),'bins',[0 1]);
%     ssfrAvgC(3)=mk_meanMedian_bin(radPosition(cmask13),(ssfr(cmask13)),'bins',[0 1]);
%     ssfrAvgC(4)=mk_ssfrMassProfs.rposMed=xmed;meanMedian_bin(radPosition(cmask14),(ssfr(cmask14)),'bins',[0 1]);
%     
%     starMassC(1)=mk_meanMedian_bin(radPosition(cmask11),(massAllGals(cmask11)),'bins',[0 1]);
%     starMassC(2)=mk_meanMedian_bin(radPosition(cmask12),(massAllGals(cmask12)),'bins',[0 1]);
%     starMassC(3)=mk_meanMedian_bin(radPosition(cmask13),(massAllGals(cmask13)),'bins',[0 1]);
%     starMassC(4)=mk_meanMedian_bin(radPosition(cmask14),(massAllGals(cmask14)),'bins',[0 1]);
%     
    

    %% save to structure
    
    
    profStruct(k).ssfrAvg=ssfrAvg;
    %profStruct(k).ssfrAvgC=ssfrAvgC;
    profStruct(k).starMass=starMass;
    %profStruct(k).starMassC=starMassC;
    
    
    
%     profStruct(k).maskStruct.mask11=mask11;
%     profStruct(k).maskStruct.mask12=mask12;
%     profStruct(k).maskStruct.mask13=mask13;
%     profStruct(k).maskStruct.mask14=mask14;
%     
%     profStruct(k).maskStruct.cmask11=cmask11;
%     profStruct(k).maskStruct.cmask12=cmask12;
%     profStruct(k).maskStruct.cmask13=cmask13;
%     profStruct(k).maskStruct.cmask14=cmask14;
end

%% find median profile over all vantage points 
quants=[0.1 0.25 0.5 0.75 0.9];
xmed=zeros(length(binEdges)-1,4);
xavg=zeros(length(binEdges)-1,4);
smassMedMed=zeros(length(quants),length(binEdges)-1,4);
smassAvgMed=zeros(length(quants),length(binEdges)-1,4);

ssfrMedMed=zeros(length(quants),length(binEdges)-1,4);
ssfrAvgMed=zeros(length(quants),length(binEdges)-1,4);


for j=1:4
    xmed(:,j)=median(squeeze(xMedian(:,:,j)),1);
    xavg(:,j)=median(squeeze(xMean(:,:,j)),1);
    
    smassMedMed(:,:,j)=quantile(squeeze(smassMedian(:,:,j)),quants,1);
    smassAvgMed(:,:,j)=quantile(squeeze(smassMean(:,:,j)),quants,1);
    
    ssfrMedMed(:,:,j)=quantile(squeeze(ssfrMedian(:,:,j)),quants,1);
    ssfrAvgMed(:,:,j)=quantile(squeeze(ssfrMean(:,:,j)),quants,1);
    
end

ssfrMassProfs.rposMed=xmed;
ssfrMassProfs.rposAvg=xavg;


ssfrMassProfs.smassMedMed=smassMedMed;
ssfrMassProfs.smassAvgMed=smassAvgMed;

ssfrMassProfs.ssfrMedMed=ssfrMedMed;
ssfrMassProfs.ssfrAvgMed=ssfrAvgMed;

ssfrMassProfs.profStruct=profStruct;
ssfrMassProfs.quants=quants;


% fprintf(['writing to: ' fname '\n']);
% save(fname,'ssfrMassProfs');    % option to save 

end
