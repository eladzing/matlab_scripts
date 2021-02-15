illustris.utils.set_illUnits(snap);

global illUnits
if readFlag
    
    fprintf('Reading in data \n');
    
    % read in data
    loadFofSub
    subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    %load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/gasProperties_snp99_TNG300.mat')
    
    % define stellar mass
    massAllGals= double(subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,:).*illUnits.massUnit); % stellar mass within 2*rhalf
    
    
    
    % set mask
    %massThresh=10^9;
    satMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'sats');
    centralMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap,'centrals');
    
    %define host mass and radius
    r200c=double(fofs.Group_R_Crit200(subsInfo.hostFof+1));
    m200c=log10(double(fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit));
    
end
%%  load obs sample 
 
% calculate distances between sub and host
global LBox
projLab={'xy' 'xz' 'yz'};
for k=1:3
    switch k
        case (1)  %% xy
            axi=[1 2];
        case (2)  %% xz
            axi=[1 3];
        case (3)  %% yz
            axi=[2 3];
    end
    
    rr(k,:)=findDistance(double(subs.SubhaloPos(axi,:)), double(fofs.GroupPos(axi,subsInfo.hostFof+1)),LBox);
    
end



radPosition=rr./r200c;
clear rr

% calculate sSFR
ssfr=double(illustris.utils.calc_ssfr(subs,'base',0));


% bin the distances
binEdges=0.1:0.2:2;
%binInd=discretize(radPosition,binEdges);
%rbin=binEdges(1:end-1)+0.5.*diff(binEdges);
%% run with all galaxies - includes ssfr=0

%bin by host mass;
mask11=satMask & m200c>=11 & m200c<12;
mask12=satMask & m200c>=12 & m200c<13;
mask13=satMask & m200c>=13 & m200c<14;
mask14=satMask & m200c>=14 & m200c<15;

% central masks
cmask11=centralMask & m200c>=11 & m200c<12;
cmask12=centralMask & m200c>=12 & m200c<13;
cmask13=centralMask & m200c>=13 & m200c<14;
cmask14=centralMask & m200c>=14 & m200c<15;



for k=1:3
    
    
    ssfrAvg(1)=mk_meanMedian_bin(radPosition(k,mask11),(ssfr(mask11)),'bins',binEdges);
    ssfrAvg(2)=mk_meanMedian_bin(radPosition(k,mask12),(ssfr(mask12)),'bins',binEdges);
    ssfrAvg(3)=mk_meanMedian_bin(radPosition(k,mask13),(ssfr(mask13)),'bins',binEdges);
    ssfrAvg(4)=mk_meanMedian_bin(radPosition(k,mask14),(ssfr(mask14)),'bins',binEdges);
    
    starMass(1)=mk_meanMedian_bin(radPosition(k,mask11),(massAllGals(mask11)),'bins',binEdges);
    starMass(2)=mk_meanMedian_bin(radPosition(k,mask12),(massAllGals(mask12)),'bins',binEdges);
    starMass(3)=mk_meanMedian_bin(radPosition(k,mask13),(massAllGals(mask13)),'bins',binEdges);
    starMass(4)=mk_meanMedian_bin(radPosition(k,mask14),(massAllGals(mask14)),'bins',binEdges);
    
    % centrals
    ssfrAvgC(1)=mk_meanMedian_bin(radPosition(k,cmask11),(ssfr(cmask11)),'bins',[0 1]);
    ssfrAvgC(2)=mk_meanMedian_bin(radPosition(k,cmask12),(ssfr(cmask12)),'bins',[0 1]);
    ssfrAvgC(3)=mk_meanMedian_bin(radPosition(k,cmask13),(ssfr(cmask13)),'bins',[0 1]);
    ssfrAvgC(4)=mk_meanMedian_bin(radPosition(k,cmask14),(ssfr(cmask14)),'bins',[0 1]);
    
    starMassC(1)=mk_meanMedian_bin(radPosition(k,cmask11),(massAllGals(cmask11)),'bins',[0 1]);
    starMassC(2)=mk_meanMedian_bin(radPosition(k,cmask12),(massAllGals(cmask12)),'bins',[0 1]);
    starMassC(3)=mk_meanMedian_bin(radPosition(k,cmask13),(massAllGals(cmask13)),'bins',[0 1]);
    starMassC(4)=mk_meanMedian_bin(radPosition(k,cmask14),(massAllGals(cmask14)),'bins',[0 1]);
    
    
    %% save to file
    
    
    profStruct(k).ssfrAvg=ssfrAvg;
    profStruct(k).ssfrAvgC=ssfrAvgC;
    profStruct(k).starMass=starMass;
    profStruct(k).starMassC=starMassC;
    profStruct(k).projection=projLab{k};
    
    
    profStruct(k).maskStruct.mask11=mask11;
    profStruct(k).maskStruct.mask12=mask12;
    profStruct(k).maskStruct.mask13=mask13;
    profStruct(k).maskStruct.mask14=mask14;
    
    profStruct(k).maskStruct.cmask11=cmask11;
    profStruct(k).maskStruct.cmask12=cmask12;
    profStruct(k).maskStruct.cmask13=cmask13;
    profStruct(k).maskStruct.cmask14=cmask14;
end
global DEFAULT_MATFILE_DIR
global simDisplayName
fname=[DEFAULT_MATFILE_DIR '/ssfr0_stellarMass_radProfiles_projected_' simDisplayName];
fprintf(['writing to: ' fname '\n']);
save(fname,'profStruct');

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
