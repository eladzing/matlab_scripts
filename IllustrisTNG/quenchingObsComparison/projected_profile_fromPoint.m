illustris.utils.set_illUnits(snap);

global illUnits
global DEFAULT_MATFILE_DIR
global simDisplayName
if readFlag
    
    fprintf('Reading in data \n');
    
    % read in data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/gasProperties_snp99_TNG300.mat')
    
    % define stellar mass
    massAllGals= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf
    
    
    
    % set mask
    %massThresh=10^9;
    %satMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'sats');
%     centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
%     
%     %define host mass and radius
     r200=double(fofs.Group_R_Mean200(subsInfo.hostFof+1)).*illUnits.lengthUnit;
%     m200c=log10(double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit));
    
end
load([DEFAULT_MATFILE_DIR '/yangSatSample_Vpoint_sig1_mean_' simDisplayName '.mat'])
fname=[DEFAULT_MATFILE_DIR '/ssfr0_stellarMass_radProfiles_projected_yangSamplePoint_' simDisplayName];

%%  load obs sample 
illustris.utils.set_illUnits(snap);

%% calculate sSFR
ssfrBase=double(illustris.utils.calc_ssfr(subs,'base',0));


% bin the distances
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

xMedian=zeros(length(satStructY),length(binEdges)-1,4);
xMean=zeros(length(satStructY),length(binEdges)-1,4);

smassMedian=zeros(length(satStructY),length(binEdges)-1,4);
smassMean=zeros(length(satStructY),length(binEdges)-1,4);

ssfrMedian=zeros(length(satStructY),length(binEdges)-1,4);
ssfrMean=zeros(length(satStructY),length(binEdges)-1,4);

for k=1:length(satStructY)
    satInd=satStructY(k).satID+1;
    hostInd=satStructY(k).hostID+1;  
    ssfr=ssfrBase(satInd);
    mass=massAllGals(satInd);
    radPosition=satStructY(k).distProj./r200c(satInd);
    hostMass=log10(double(fofs.Group_M_Crit200(hostInd)).*illUnits.massUnit);
    
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

fprintf(['writing to: ' fname '\n']);
save(fname,'ssfrMassProfs');

% %% add in the obs data
% global DEFAULT_MATFILE_DIR
% load([DEFAULT_MATFILE_DIR '/ssfr_rpos_dataGrab.mat'])
%
%
% %% plot
% global simDisplayName
% cc=brewermap(8,'Set1');
% h=[];
% hf=myFigure; % figure('position',[953   261   938   737])
%
% for i=1:4
%     yy=log10(ssfrAvg(i,:,3));
%
%     switch i
%         case 1
%             colo=cc(2,:);
%             nam='11-12, TNG';
%         case 2
%             colo=cc(5,:);
%             nam='12-13, TNG';
%
%         case 3
%             colo=cc(3,:);
%             nam='13-14, TNG';
%
%         case 4
%             colo=cc(1,:);
%             nam='14-15, TNG';
%
%     end
%     h(i)=plot(rbin,yy,'color',colo,'DisplayName',nam,"LineStyle","--");
%     if i==1; hold on; end
%     plot(0,log10(ssfrAvgC(i)),'s','color',colo)
% end
%
% % obs
% for i=5:8
%
%     switch i
%         case 5
%             xx=ssfr_rpos_11_main(2,:);
%             yy=ssfr_rpos_11_main(1,:);
%             pos=ssfr_rpos_11_top(1,:)-ssfr_rpos_11_main(1,:);
%             neg=ssfr_rpos_11_main(1,:)-ssfr_rpos_11_bottom(1,:);
%             colo=cc(2,:);
%             nam='11-12, Obs';
%
%             ccc=ssfr_rpos_11_central(1,2);
%             cp=ssfr_rpos_11_central(1,3)-ssfr_rpos_11_central(1,2);
%             cm=ssfr_rpos_11_central(1,2)-ssfr_rpos_11_central(1,1);
%         case 6
%             xx=ssfr_rpos_12_main(2,:);
%             yy=ssfr_rpos_12_main(1,:);
%             pos=ssfr_rpos_12_top(1,:)-ssfr_rpos_12_main(1,:);
%             neg=ssfr_rpos_12_main(1,:)-ssfr_rpos_12_bottom(1,:);
%             colo=cc(5,:);
%             nam='12-13, Obs';
%
%             ccc=ssfr_rpos_12_central(1,2);
%             cp=ssfr_rpos_12_central(1,3)-ssfr_rpos_12_central(1,2);
%             cm=ssfr_rpos_12_central(1,2)-ssfr_rpos_12_central(1,1);
%
%         case 7
%             xx=ssfr_rpos_13_main(2,:);
%             yy=ssfr_rpos_13_main(1,:);
%             pos=ssfr_rpos_13_top(1,:)-ssfr_rpos_13_main(1,:);
%             neg=ssfr_rpos_13_main(1,:)-ssfr_rpos_13_bottom(1,:);
%             colo=cc(3,:);
%             nam='13-14, Obs';
%
%             ccc=ssfr_rpos_13_central(1,2);
%             cp=ssfr_rpos_13_central(1,3)-ssfr_rpos_13_central(1,2);
%             cm=ssfr_rpos_13_central(1,2)-ssfr_rpos_13_central(1,1);
%
%         case 8
%             xx=ssfr_rpos_14_main(2,:);
%             yy=ssfr_rpos_14_main(1,:);
%             pos=ssfr_rpos_14_top(1,:)-ssfr_rpos_14_main(1,:);
%             neg=ssfr_rpos_14_main(1,:)-ssfr_rpos_14_bottom(1,:);
%             colo=cc(1,:);
%             nam='14-15, Obs';
%
%             ccc=ssfr_rpos_14_central(1,2);
%             cp=ssfr_rpos_14_central(1,3)-ssfr_rpos_14_central(1,2);
%             cm=ssfr_rpos_14_central(1,2)-ssfr_rpos_14_central(1,1);
%
%     end
%
%
%
%
%     h(i)=errorbar(xx,yy,neg,pos,'color',colo,...
%         'DisplayName',nam);
%
%     errorbar(0,ccc,cm,cp,'x','color',colo)
% end
%
% set(gca,'fontsize',14)
% %hl=legend(h,'location','northwest' );
% hl=legend(h,'numcolumns',2,'location','northwest' );
%
%
% xlim([-0.06 1.5])
% ylim([-12.05 -9])
% xlabelmine('$r/R_\mathrm{200,c}$');
% ylabelmine('sSFR');
% titlemine(['Comparison with ' simDisplayName]);
%
%
% %printout_fig(gcf,['obs_compare_ssfr_profile_' simDisplayName])
%
%
