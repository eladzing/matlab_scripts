function res= generate_hih2_profilesFull_3d(hih2Struct,fofs,subs,varargin)

global illUnits
% global DEFAULT_MATFILE_DIR
% global simDisplayName


virType='crit200';
sampleMask=hih2Struct.galMask;
radBinEdges=0.1:0.1:2;
sfrBinEdges=-13:0.5:-8;
mainSequenceType='all';
modelType='br';
hiMassThresh=100; % lower enforced limit to avoid NAN's in log for zero mass
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
        case {'radbinedges','radbins'}
            i=i+1;
            radBinEdges=varargin{i};
        case {'sfrbins','sfrbinedges'}
            i=i+1;
            sfrBinEdges=varargin{i};
        case 'mask'
            i=i+1;
            sampleMask=varargin{i};
        case 'mainsequence'
            i=i+1;
            mainSequenceType=varargin{i};
        case 'model'
            i=i+1;
            modelType=varargin{i};
        case {'massthresh', 'massthreshold'}
            i=i+1;
            hiMassThresh=varargin{i};
        case {'br','gk','kmt'}
            modelType=varargin{i};
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end

switch lower(modelType)
    case 'br'
        kmodel=1;
    case 'gk'
        kmodel=2;
    case 'kmt'
        kmodel=3;
end

%% Perliminaries

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
% define stellar mass
gMass= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf

%% generate 'main-sequence'
switch lower(mainSequenceType)
    case 'central'
        mmaskG=hih2Struct.galMask & hih2Struct.Gal.GalHIMass(kmodel,:)>0 & subsInfo.isCentral;
        mmaskC=hih2Struct.galMask & hih2Struct.CGMall.CGMallHIMass(kmodel,:)>0 & subsInfo.isCentral;
    case 'all'
        mmaskG=hih2Struct.galMask & hih2Struct.Gal.GalHIMass(kmodel,:)>0;
        mmaskC=hih2Struct.galMask & hih2Struct.CGMall.CGMallHIMass(kmodel,:)>0;
end

galH_star=mk_meanMedian_bin(log10(hih2Struct.galMass(mmaskG)),...
    log10(hih2Struct.Gal.GalHIMass(kmodel,mmaskG)./hih2Struct.galMass(mmaskG)),'bins',8.29:0.05:12.5);


cgmH_star=mk_meanMedian_bin(log10(hih2Struct.galMass(mmaskC)),...
    log10(hih2Struct.CGMall.CGMallHIMass(kmodel,mmaskC)./hih2Struct.galMass(mmaskC)),'bins',8.29:0.05:12.5);


%% calculate sSFR
ssfrBase=double(illustris.utils.calc_ssfr(subs,'base',1e-13));
ssfr=ssfrBase(sampleMask);

%% set mask
loc.Gal=mask_structure(hih2Struct.Gal,sampleMask);
loc.CGMin=mask_structure(hih2Struct.CGMin,sampleMask);
loc.CGMout=mask_structure(hih2Struct.CGMout,sampleMask);
loc.CGMall=mask_structure(hih2Struct.CGMall,sampleMask);
loc.Sub=mask_structure(hih2Struct.Sub,sampleMask);
comps=fields(loc);
gMass=gMass(sampleMask);

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
hmask(1,:)=log10(mvir)>=10.7 & log10(mvir)<12;
hmask(2,:)=log10(mvir)>=12 & log10(mvir)<13;
hmask(3,:)=log10(mvir)>=13 & log10(mvir)<14;
hmask(4,:)=log10(mvir)>=14 & log10(mvir)<15;

res.hostNums=sum(hmask,2);
res.hostHeader='Host Mass bins: 10.7, 12, 13, 14, 15';

%% split by stellarmass - 4 groups
lgMass=log10(gMass);
smask4=false(4,length(gMass));
smask4(1,:)=lgMass>=8 & lgMass<9;
smask4(2,:)=lgMass>=9 & lgMass<10;
smask4(3,:)=lgMass>=10 & lgMass<11;
smask4(4,:)=lgMass>=11 ;
res.byHostStar4.massHeader='stellar mass bins: 8 9 10 11 and above';
%% split by stellarmass - 2 groups
smask2=false(2,length(gMass));
smask2(1,:)=lgMass>=8 & lgMass<10.5;
smask2(2,:)=lgMass>=10.5;
res.byHostStar2.massHeader='stellar mass bins: 8 10.5 and above';

%% by Host mass only
for j=1:4
    mask=squeeze(hmask(j,:));
    ssfrProf=mk_meanMedian_bin(radPosition(mask),ssfr(mask),'bins',radBinEdges);
    starMassProf=mk_meanMedian_bin(radPosition(mask),gMass(mask),'bins',radBinEdges);
    
    byHost.ssfrAvg(:,j)=ssfrProf.yMean;
    byHost.ssfrMed(:,j)=ssfrProf.yMedian;
    
    byHost.gMassAvg(:,j)=starMassProf.yMean;
    byHost.gMassMed(:,j)=starMassProf.yMedian;
    
    byHost.xMed(:,j)=starMassProf.xMedian;
    byHost.xAvg(:,j)=starMassProf.xMean;
    byHost.count(:,j)=starMassProf.binCount;
    
    for k=1:length(comps)
        hiMass=loc.(comps{k}).(strcat(comps{k},'HIMass'));
        h2Mass=loc.(comps{k}).(strcat(comps{k},'H2Mass'));
        for jj=1:3
            hiProf=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask),'bins',radBinEdges);
            h2Prof=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask),'bins',radBinEdges);
            
            byHost.(comps{k}).hiMassAvg(jj,:,j)=hiProf.yMean;
            byHost.(comps{k}).hiMassMed(jj,:,j)=hiProf.yMedian;
            
            byHost.(comps{k}).h2MassAvg(jj,:,j)=h2Prof.yMean;
            byHost.(comps{k}).h2MassMed(jj,:,j)=h2Prof.yMedian;
            
            % normalized by stellar mass
            hiProfN=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask)./gMass(mask),'bins',radBinEdges);
            h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',radBinEdges);
            
            byHost.(comps{k}).hiMassAvgN(jj,:,j)=hiProfN.yMean;
            byHost.(comps{k}).hiMassMedN(jj,:,j)=hiProfN.yMedian;
            
            byHost.(comps{k}).h2MassAvgN(jj,:,j)=h2ProfN.yMean;
            byHost.(comps{k}).h2MassMedN(jj,:,j)=h2ProfN.yMedian;
            
            %% hi offset
            if strcmp(comps{k},'Gal')
                ms=galH_star.xMedian;
                hm=galH_star.yMedian;
            elseif strcmp(comps{k},'CGMall')
                ms=cgmH_star.xMedian;
                hm=cgmH_star.yMedian;
            end
            ms=ms(~isnan(ms));
            hm=hm(~isnan(ms));
            msMass=interp1(ms,hm,log10(gMass(mask)));
            
            him=hiMass(jj,mask);
            him(him<hiMassThresh)=hiMassThresh;
            hiDef=-1.*(him./gMass(mask)-10.^msMass);
            hiDefL=-1.*(log10(him./gMass(mask))-msMass);
            
            hiProfD=mk_meanMedian_bin(radPosition(mask),hiDef,'bins',radBinEdges);
            hiProfDL=mk_meanMedian_bin(radPosition(mask),hiDefL,'bins',radBinEdges);
            %hiProfDN=mk_meanMedian_bin(radPosition(mask),-1.*(him./gMass(mask))./(10.^msMass)-1,'bins',radBinEdges);
            
            %h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',binEdges);
            
            byHost.(comps{k}).hiMassAvgD(jj,:,j)=hiProfD.yMean;
            byHost.(comps{k}).hiMassMedD(jj,:,j)=hiProfD.yMedian;
            
            byHost.(comps{k}).hiMassAvgDL(jj,:,j)=hiProfDL.yMean;
            byHost.(comps{k}).hiMassMedDL(jj,:,j)=hiProfDL.yMedian;
            
            %byHost.(comps{k}).h2MassAvgN(jj,:,j)=h2ProfN.yMean;
            %byHost.(comps{k}).h2MassMedN(jj,:,j)=h2ProfN.yMedian;
            
            %% mass  by ssfr
            hiProfSFRN=mk_meanMedian_bin(log10(ssfr(mask)),hiMass(jj,mask)./gMass(mask),'bins',sfrBinEdges);
            h2ProfSFRN=mk_meanMedian_bin(log10(ssfr(mask)),h2Mass(jj,mask)./gMass(mask),'bins',sfrBinEdges);
            
            byHost.xsfrMed(:,j)=hiProfSFRN.xMedian;
            byHost.xsfrAvg(:,j)=hiProfSFRN.xMean;
            
            byHost.(comps{k}).hiMassSFRAvgN(jj,:,j)=hiProfSFRN.yMean;
            byHost.(comps{k}).hiMassSFRMedN(jj,:,j)=hiProfSFRN.yMedian;
            
            byHost.(comps{k}).h2MassSFRAvgN(jj,:,j)=h2ProfSFRN.yMean;
            byHost.(comps{k}).h2MassSFRMedN(jj,:,j)=h2ProfSFRN.yMedian;
            
            
            %% deficiency by ssfr
            hiProfDSFR=mk_meanMedian_bin(log10(ssfr(mask)),hiDef,'bins',sfrBinEdges);
            hiProfDLSFR=mk_meanMedian_bin(log10(ssfr(mask)),hiDefL,'bins',sfrBinEdges);
            
            %byHost.xsfrMed(:,j)=hiProfSFRN.xMedian;
            %byHost.xsfrAvg(:,j)=hiProfSFRN.xMean;
            
            byHost.(comps{k}).hiDefSFRAvg(jj,:,j)=hiProfDSFR.yMean;
            byHost.(comps{k}).hiDefSFRMed(jj,:,j)=hiProfDSFR.yMedian;
            
            byHost.(comps{k}).hiDefLSFRAvg(jj,:,j)=hiProfDLSFR.yMean;
            byHost.(comps{k}).hiDefLSFRMed(jj,:,j)=hiProfDLSFR.yMedian;
            
        end
    end
end

%% by host and stellar mask (2 groups)

for i=1:2
    
    for j=1:4
        mask=squeeze(hmask(j,:)) & squeeze(smask2(i,:));
        ssfrProf=mk_meanMedian_bin(radPosition(mask),ssfr(mask),'bins',radBinEdges);
        starMassProf=mk_meanMedian_bin(radPosition(mask),gMass(mask),'bins',radBinEdges);
        
        byHostStar2(i).ssfrAvg(:,j)=ssfrProf.yMean;
        byHostStar2(i).ssfrMed(:,j)=ssfrProf.yMedian;
        
        byHostStar2(i).gMassAvg(:,j)=starMassProf.yMean;
        byHostStar2(i).gMassMed(:,j)=starMassProf.yMedian;
        
        byHostStar2(i).xMed(:,j)=starMassProf.xMedian;
        byHostStar2(i).xAvg(:,j)=starMassProf.xMean;
        
        byHostStar2(i).count(:,j)=starMassProf.binCount;
        
        
        for k=1:length(comps)
            hiMass=loc.(comps{k}).(strcat(comps{k},'HIMass'));
            h2Mass=loc.(comps{k}).(strcat(comps{k},'H2Mass'));
            for jj=1:3
                hiProf=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask),'bins',radBinEdges);
                h2Prof=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask),'bins',radBinEdges);
                
                byHostStar2(i).(comps{k}).hiMassAvg(jj,:,j)=hiProf.yMean;
                byHostStar2(i).(comps{k}).hiMassMed(jj,:,j)=hiProf.yMedian;
                
                byHostStar2(i).(comps{k}).h2MassAvg(jj,:,j)=h2Prof.yMean;
                byHostStar2(i).(comps{k}).h2MassMed(jj,:,j)=h2Prof.yMedian;
                
                % normalized by stellar mass
                hiProfN=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask)./gMass(mask),'bins',radBinEdges);
                h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',radBinEdges);
                
                byHostStar2(i).(comps{k}).hiMassAvgN(jj,:,j)=hiProfN.yMean;
                byHostStar2(i).(comps{k}).hiMassMedN(jj,:,j)=hiProfN.yMedian;
                
                byHostStar2(i).(comps{k}).h2MassAvgN(jj,:,j)=h2ProfN.yMean;
                byHostStar2(i).(comps{k}).h2MassMedN(jj,:,j)=h2ProfN.yMedian;
                
                %%hi offset
                if strcmp(comps{k},'Gal')
                    ms=galH_star.xMedian;
                    hm=galH_star.yMedian;
                elseif strcmp(comps{k},'CGMall')
                    ms=cgmH_star.xMedian;
                    hm=cgmH_star.yMedian;
                end
                ms=ms(~isnan(ms));
                hm=hm(~isnan(ms));
                msMass=interp1(ms,hm,log10(gMass(mask)));
                
                him=hiMass(jj,mask);
                him(him<hiMassThresh)=hiMassThresh;
                hiDef=-1.*(him./gMass(mask)-10.^msMass);
                hiDefL=-1.*(log10(him./gMass(mask))-msMass);
                
                hiProfD=mk_meanMedian_bin(radPosition(mask),hiDef,'bins',radBinEdges);
                hiProfDL=mk_meanMedian_bin(radPosition(mask),hiDefL,'bins',radBinEdges);
                %hiProfDN=mk_meanMedian_bin(radPosition(mask),-1.*(him./gMass(mask))./(10.^msMass)-1,'bins',radBinEdges);
                
                %h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',binEdges);
                
                byHostStar2(i).(comps{k}).hiMassAvgD(jj,:,j)=hiProfD.yMean;
                byHostStar2(i).(comps{k}).hiMassMedD(jj,:,j)=hiProfD.yMedian;
                
                byHostStar2(i).(comps{k}).hiMassAvgDL(jj,:,j)=hiProfDL.yMean;
                byHostStar2(i).(comps{k}).hiMassMedDL(jj,:,j)=hiProfDL.yMedian;
                
                %byHostStar2(i).(comps{k}).h2MassAvgN(jj,:,j)=h2ProfN.yMean;
                %byHostStar2(i).(comps{k}).h2MassMedN(jj,:,j)=h2ProfN.yMedian;
                
                %% mass by ssfr
                hiProfSFRN=mk_meanMedian_bin(log10(ssfr(mask)),hiMass(jj,mask)./gMass(mask),'bins',sfrBinEdges);
                h2ProfSFRN=mk_meanMedian_bin(log10(ssfr(mask)),h2Mass(jj,mask)./gMass(mask),'bins',sfrBinEdges);
                
                byHostStar2(i).xsfrMed(:,j)=hiProfSFRN.xMedian;
                byHostStar2(i).xsfrAvg(:,j)=hiProfSFRN.xMean;
                
                byHostStar2(i).(comps{k}).hiMassSFRAvgN(jj,:,j)=hiProfSFRN.yMean;
                byHostStar2(i).(comps{k}).hiMassSFRMedN(jj,:,j)=hiProfSFRN.yMedian;
                
                byHostStar2(i).(comps{k}).h2MassSFRAvgN(jj,:,j)=h2ProfSFRN.yMean;
                byHostStar2(i).(comps{k}).h2MassSFRMedN(jj,:,j)=h2ProfSFRN.yMedian;
                
                
                %% deficiency by ssfr
                hiProfDSFR=mk_meanMedian_bin(log10(ssfr(mask)),hiDef,'bins',sfrBinEdges);
                hiProfDLSFR=mk_meanMedian_bin(log10(ssfr(mask)),hiDefL,'bins',sfrBinEdges);
                
                byHostStar2(i).(comps{k}).hiDefSFRAvg(jj,:,j)=hiProfDSFR.yMean;
                byHostStar2(i).(comps{k}).hiDefSFRMed(jj,:,j)=hiProfDSFR.yMedian;
                
                byHostStar2(i).(comps{k}).hiDefLSFRAvg(jj,:,j)=hiProfDLSFR.yMean;
                byHostStar2(i).(comps{k}).hiDefLSFRMed(jj,:,j)=hiProfDLSFR.yMedian;
                
                
            end
        end
    end
end


%% by host and stellar mass 4 groups

for i=1:4  % in stellar mass bins
    
    for j=1:4 % in host mass bins
        mask=squeeze(hmask(j,:)) & squeeze(smask4(i,:));
        ssfrProf=mk_meanMedian_bin(radPosition(mask),ssfr(mask),'bins',radBinEdges);
        starMassProf=mk_meanMedian_bin(radPosition(mask),gMass(mask),'bins',radBinEdges);
        
        byHostStar4(i).ssfrAvg(:,j)=ssfrProf.yMean;
        byHostStar4(i).ssfrMed(:,j)=ssfrProf.yMedian;
        
        byHostStar4(i).gMassAvg(:,j)=starMassProf.yMean;
        byHostStar4(i).gMassMed(:,j)=starMassProf.yMedian;
        
        byHostStar4(i).xMed(:,j)=starMassProf.xMedian;
        byHostStar4(i).xAvg(:,j)=starMassProf.xMean;
        
        byHostStar4(i).count(:,j)=starMassProf.binCount;
        
        for k=1:length(comps)
            hiMass=loc.(comps{k}).(strcat(comps{k},'HIMass'));
            h2Mass=loc.(comps{k}).(strcat(comps{k},'H2Mass'));
            for jj=1:3
                hiProf=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask),'bins',radBinEdges);
                h2Prof=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask),'bins',radBinEdges);
                
                byHostStar4(i).(comps{k}).hiProf(jj,j)=hiProf;
                byHostStar4(i).(comps{k}).h2Prof(jj,j)=h2Prof;
                
                
                % normalized by stellar mass
                hiProfN=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask)./gMass(mask),'bins',radBinEdges);
                h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',radBinEdges);

                byHostStar4(i).(comps{k}).hiProfN(jj,j)=hiProfN;
                byHostStar4(i).(comps{k}).h2ProfN(jj,j)=h2ProfN;

                               
                
                %hi offset
                if strcmp(comps{k},'Gal')
                    ms=galH_star.xMedian;
                    hm=galH_star.yMedian;
                elseif strcmp(comps{k},'CGMall')
                    ms=cgmH_star.xMedian;
                    hm=cgmH_star.yMedian;
                end
                ms=ms(~isnan(ms));
                hm=hm(~isnan(ms));
                msMass=interp1(ms,hm,log10(gMass(mask)));
                
                him=hiMass(jj,mask);
                him(him<hiMassThresh)=hiMassThresh;
                hiDef=-1.*(him./gMass(mask)-10.^msMass);
                hiDefL=-1.*(log10(him./gMass(mask))-msMass);
                
                hiProfD=mk_meanMedian_bin(radPosition(mask),hiDef,'bins',radBinEdges);
                hiProfDL=mk_meanMedian_bin(radPosition(mask),hiDefL,'bins',radBinEdges);
                %hiProfDN=mk_meanMedian_bin(radPosition(mask),-1.*(him./gMass(mask))./(10.^msMass)-1,'bins',radBinEdges);
                
                %h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',binEdges);
                
                byHostStar4(i).(comps{k}).hiProfD(jj,j)=hiProfD;
                byHostStar4(i).(comps{k}).hiProfDL(jj,j)=hiProfDL;
                
                
                %% mass by ssfr
                hiProfSFRN=mk_meanMedian_bin(log10(ssfr(mask)),hiMass(jj,mask)./gMass(mask),'bins',sfrBinEdges);
                h2ProfSFRN=mk_meanMedian_bin(log10(ssfr(mask)),h2Mass(jj,mask)./gMass(mask),'bins',sfrBinEdges);
                
                byHostStar4(i).(comps{k}).hiProfSFRN(jj,j)=hiProfSFRN;
                byHostStar4(i).(comps{k}).h2ProfSFRN(jj,j)=h2ProfSFRN;

                %% deficiency by ssfr
                hiProfDSFR=mk_meanMedian_bin(log10(ssfr(mask)),hiDef,'bins',sfrBinEdges);
                hiProfDLSFR=mk_meanMedian_bin(log10(ssfr(mask)),hiDefL,'bins',sfrBinEdges);
                
                byHostStar4(i).(comps{k}).hiProfDSFRN(jj,j)=hiProfDSFR;
                byHostStar4(i).(comps{k}).hiProfDLSFRN(jj,j)=hiProfDLSFR;
                                
                
            end
        end
    end
end


















% %% split by stellarmass
% smask=false(3,length(gMass));
% smask(1,:)=log10(gMass)>=9 & log10(gMass)<10;
% smask(2,:)=log10(gMass)>=10 & log10(gMass)<11;
% smask(3,:)=log10(gMass)>=11 ;
%
% for j=1:3
%     mask=squeeze(smask(j,:));
%     ssfrProf=mk_meanMedian_bin(radPosition(mask),ssfr(mask),'bins',binEdges);
%     starMassProf=mk_meanMedian_bin(radPosition(mask),gMass(mask),'bins',binEdges);
%
%     byMass.ssfrAvg(:,j)=ssfrProf.yMean;
%     byMass.ssfrMed(:,j)=ssfrProf.yMedian;
%
%     byMass.gMassAvg(:,j)=starMassProf.yMean;
%     byMass.gMassMed(:,j)=starMassProf.yMedian;
%
%     byMass.xMed(:,j)=starMassProf.xMedian;
%     byMass.xAvg(:,j)=starMassProf.xMean;
%
%     for k=1:length(comps)
%         hiMass=loc.(comps{k}).(strcat(comps{k},'HIMass'));
%         h2Mass=loc.(comps{k}).(strcat(comps{k},'H2Mass'));
%         for jj=1:3
%             hiProf=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask),'bins',binEdges);
%             h2Prof=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask),'bins',binEdges);
%
%             byMass.(comps{k}).hiMassAvg(jj,:,j)=hiProf.yMean;
%             byMass.(comps{k}).hiMassMed(jj,:,j)=hiProf.yMedian;
%
%             byMass.(comps{k}).h2MassAvg(jj,:,j)=h2Prof.yMean;
%             byMass.(comps{k}).h2MassMed(jj,:,j)=h2Prof.yMedian;
%
%             % normalized by stellar mass
%             hiProfN=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask)./gMass(mask),'bins',binEdges);
%             h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',binEdges);
%
%             byMass.(comps{k}).hiMassAvgN(jj,:,j)=hiProfN.yMean;
%             byMass.(comps{k}).hiMassMedN(jj,:,j)=hiProfN.yMedian;
%
%             byMass.(comps{k}).h2MassAvgN(jj,:,j)=h2ProfN.yMean;
%             byMass.(comps{k}).h2MassMedN(jj,:,j)=h2ProfN.yMedian;
%
%
% %              %%hi offset
% %             if strcmp(comps{k},'Gal')
% %                 ms=galH_star.xMedian;
% %                 hm=galH_star.yMedian;
% %             elseif strcmp(comps{k},'CGMall')
% %                 ms=cgmH_star.xMedian;
% %                 hm=cgmH_star.yMedian;
% %             end
% %             ms=ms(~isnan(ms));
% %             hm=hm(~isnan(ms));
% %             msMass=interp1(ms,hm,log10(gMass(mask)));
% %
% %             hiProfD=mk_meanMedian_bin(radPosition(mask),(hiMass(jj,mask)-10.^msMass),'bins',binEdges);
% %             hiProfDN=mk_meanMedian_bin(radPosition(mask),hiMass(jj,mask)./(10.^msMass)-1,'bins',binEdges);
% %             %h2ProfN=mk_meanMedian_bin(radPosition(mask),h2Mass(jj,mask)./gMass(mask),'bins',binEdges);
% %
% %             byMass.(comps{k}).hiMassAvgD(jj,:,j)=hiProfD.yMean;
% %             byMass.(comps{k}).hiMassMedD(jj,:,j)=hiProfD.yMedian;
% %
% %             byMass.(comps{k}).hiMassAvgDN(jj,:,j)=hiProfDN.yMean;
% %             byMass.(comps{k}).hiMassMedDN(jj,:,j)=hiProfDN.yMedian;
%
%         end
%     end
% end
%
%%
res.byHost=byHost;
res.byHostStar2=byHostStar2;
res.byHostStar4=byHostStar4;
%res.byMass=byMass;

end
