function res= generate_hih2_profiles_3d(hih2Struct,fofs,subs,varargin)

global illUnits
% global DEFAULT_MATFILE_DIR
% global simDisplayName


virType='crit200';
sampleMask=hih2Struct.galMask;
binEdges=0.1:0.1:2;
%% parse arguments
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'crit','crit200','200crit'}
            virType='crit200';
        case {'mean','mean200','200mean'}
            virType='mean200';
        case {'500','crit500','500crit'}
            virType='crit500';
        case {'binedges','bins'}
            i=i+1;
            binEdges=varargin{i};
        case 'mask'
            i=i+1;
            sampleMask=varargin{i};
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end

%% Perliminaries

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
% define stellar mass
gMass= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf


%% calculate sSFR
ssfrBase=double(illustris.utils.calc_ssfr(subs,'base',0));
ssfr=ssfrBase(sampleMask);
%% set mask
loc.Gal=mask_structure(hih2Struct.Gal,sampleMask);
loc.CGMin=mask_structure(hih2Struct.CGMin,sampleMask);
loc.CGMout=mask_structure(hih2Struct.CGMout,sampleMask);
loc.CGMall=mask_structure(hih2Struct.CGMall,sampleMask);
loc.Sub=mask_structure(hih2Struct.Sub,sampleMask);
comps=fields(loc);
gMass=gMass(sampleMask);

list
%% define virial parameters
switch lower(virType)
    case 'crit200'
        rvir=double(fofs.Group_R_Crit200(subsInfo.hostFof(sampleMask)+1)); % host rvir for each galaxy
        mvir=double(fofs.Group_M_Crit200(subsInfo.hostFof(sampleMask)+1)).*illUnits.massUnit;  % host mvir for each galaxy
    case 'mean200'
        rvir=double(fofs.Group_R_Mean200(subsInfo.hostFof(sampleMask)+1));
        mvir=double(fofs.Group_M_Mean200(subsInfo.hostFof(sampleMask)+1)).*illUnits.massUnit;
    case 'crit500'
        rvir=double(fofs.Group_R_Crit500(subsInfo.hostFof(sampleMask)+1));
        mvir=double(fofs.Group_M_Crit500(subsInfo.hostFof(sampleMask)+1)).*illUnits.massUnit;
end
%     m200c=log10(double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit));


%% calculate distances
global LBox
radPosition = findDistance(subs.SubhaloPos(:,sampleMask),fofs.GroupPos(:,subsInfo.hostFof(sampleMask)+1),...
    LBox,3)./rvir;

%% initialize


%% split by hostmass
hmask=false(4,length(mvir));
hmask(1,:)=log10(mvir)>=11 & log10(mvir)<12;
hmask(2,:)=log10(mvir)>=12 & log10(mvir)<13;
hmask(3,:)=log10(mvir)>=13 & log10(mvir)<14;
hmask(4,:)=log10(mvir)>=14 & log10(mvir)<15;


for j=1:4
    mask=squeeze(hmask(j,:));
    ssfrProf=mk_meanMedian_bin(radPosition(mask),ssfr(mask),'bins',binEdges);
    starMassProf=mk_meanMedian_bin(radPosition(mask),gMass(mask),'bins',binEdges);
    
    byHost.ssfrAvg(:,j)=ssfrProf.yMean;
    byHost.ssfrMed(:,j)=ssfrProf.yMedian;
    
    byHost.gMassAvg(:,j)=starMassProf.yMean;
    byHost.gMassMed(:,j)=starMassProf.yMedian;
    
    byHost.xMed(:,j)=starMassProf.xMedian;
    byHost.xAvg(:,j)=starMassProf.xMean;
    
    for k=1:length{comps}
        hiMass=loc.(comps{k}).(strcat(comps{k},'HIMass'));
        h2Mass=loc.(comps{k}).(strcat(comps{k},'H2Mass'));
        for jj=1:3
            hiProf=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask),'bins',binEdges);
            h2Prof=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask),'bins',binEdges);
            
            byHost.(comps{k}).hiMassAvg(jj,:,j)=hiProf.yMean;
            byHost.(comps{k}).hiMassMed(jj,:,j)=hiProf.yMedian;
            
            byHost.(comps{k}).h2MassAvg(jj,:,j)=h2Prof.yMean;
            byHost.(comps{k}).h2MassMed(jj,:,j)=h2Prof.yMedian;
        end
    end
end

%% split by stellarmass
smask=false(3,length(gMass));
smask(1,:)=log10(gMass)>=9 & log10(gMass)<10;
smask(2,:)=log10(gMass)>=10 & log10(gMass)<11;
smask(3,:)=log10(gMass)>=11 ;

for j=1:3
    mask=squeeze(smask(j,:));
    ssfrProf=mk_meanMedian_bin(radPosition(mask),ssfr(mask),'bins',binEdges);
    starMassProf=mk_meanMedian_bin(radPosition(mask),gMass(mask),'bins',binEdges);
    
    byMass.ssfrAvg(:,j)=ssfrProf.yMean;
    byMass.ssfrMed(:,j)=ssfrProf.yMedian;
    
    byMass.gMassAvg(:,j)=starMassProf.yMean;
    byMass.gMassMed(:,j)=starMassProf.yMedian;
    
    byMass.xMed(:,j)=starMassProf.xMedian;
    byMass.xAvg(:,j)=starMassProf.xMean;
    
    for k=1:length{comps}
        hiMass=loc.(comps{k}).(strcat(comps{k},'HIMass'));
        h2Mass=loc.(comps{k}).(strcat(comps{k},'H2Mass'));
        for jj=1:3
            hiProf=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask),'bins',binEdges);
            h2Prof=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask),'bins',binEdges);
            
            byMass.(comps{k}).hiMassAvg(jj,:,j)=hiProf.yMean;
            byMass.(comps{k}).hiMassMed(jj,:,j)=hiProf.yMedian;
            
            byMass.(comps{k}).h2MassAvg(jj,:,j)=h2Prof.yMean;
            byMass.(comps{k}).h2MassMed(jj,:,j)=h2Prof.yMedian;
        end
    end
end
    
%%
res.byHost=byHost;
res.byMass=byMass;
        
    end
