%% plot the number of gas cells (without star-forming gas) 
%  within the Gal, CGM, Out, components, also seperating into SF and Quenched galaxies 
snap=99;
global simDisplayName
global illUnits
global DEFAULT_MATFILE_DIR
if readFlag
load([DEFAULT_MATFILE_DIR '/gasProperties_snp' num2str(snap) '_' simDisplayName '.mat'])

global DRACOFLAG
        if DRACOFLAG
            fofs=illustris.groupcat.loadHalos(bp,snap);
            subs=illustris.groupcat.loadSubhalos(bp,snap);
    
        else
            loadFofSubTNG100
        end
end
subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);


galMass=tCoolStruct.galMass(tCoolStruct.galMask);  % galaxy stellar mass
ssfr=illustris.utils.calc_ssfr(subs);% sfr./galMass + 10^0.5*1e-17.*10.^(0.5.*rand(size(sfr)));
ssfr=ssfr(tCoolStruct.galMask);
%gasTemp=tCoolStruct.inGal.meanTempMW(1,tCoolStruct.galMask);

ncGalBase=tCoolStruct.inGal.cellNum(tCoolStruct.galMask);
ncCGMBase=tCoolStruct.inCGM.cellNum(tCoolStruct.galMask);
ncOutBase=tCoolStruct.inOut.cellNum(tCoolStruct.galMask);

%% make masks 
centralMask= subsInfo.isCentral(tCoolStruct.galMask);
spbMask = mk_splashback_mask('time',5,'both',0.1);

ssfrThresh=1e-11;
qMask=ssfr<=ssfrThresh;

galaxyMask=centralMask & ~spbMask(tCoolStruct.galMask);

ncGal=ncGalBase(galaxyMask);
ncCGM=ncCGMBase(galaxyMask);
ncOut=ncOutBase(galaxyMask);

cellNumStruct.inGal.ncQuant=quantile(ncGal,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inCGM.ncQuant=quantile(ncCGM,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inOut.ncQuant=quantile(ncOut,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inGal.ncAvg=mean(ncGal);
cellNumStruct.inCGM.ncAvg=mean(ncCGM);
cellNumStruct.inOut.ncAvg=mean(ncOut);

ncGalSF=ncGalBase(galaxyMask & ~qMask);
ncCGMSF=ncCGMBase(galaxyMask & ~qMask);
ncOutSF=ncOutBase(galaxyMask & ~qMask);

cellNumStruct.inGal.ncQuantSF=quantile(ncGalSF,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inCGM.ncQuantSF=quantile(ncCGMSF,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inOut.ncQuantSF=quantile(ncOutSF,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inGal.ncAvgSF=mean(ncGalSF);
cellNumStruct.inCGM.ncAvgSF=mean(ncCGMSF);
cellNumStruct.inOut.ncAvgSF=mean(ncOutSF);

ncGalQ=ncGalBase(galaxyMask & qMask);
ncCGMQ=ncCGMBase(galaxyMask & qMask);
ncOutQ=ncOutBase(galaxyMask & qMask);

cellNumStruct.inGal.ncQuantQ=quantile(ncGalQ,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inCGM.ncQuantQ=quantile(ncCGMQ,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inOut.ncQuantQ=quantile(ncOutQ,[0.1 0.25 0.5 0.75 0.9]);
cellNumStruct.inGal.ncAvgQ=mean(ncGalQ);
cellNumStruct.inCGM.ncAvgQ=mean(ncCGMQ);
cellNumStruct.inOut.ncAvgQ=mean(ncOutQ);

hostmBase=fofs.Group_M_Crit200(subsInfo.hostFof+1).*illUnits.massUnit;
hostm=hostmBase(tCoolStruct.galMask);
hostmm=hostm(galaxyMask);

ncGMH=mk_meanMedian_bin(log10(hostmm),ncGal,'nbins',20);
ncCMH=mk_meanMedian_bin(log10(hostmm),ncCGM,'nbins',20);
ncOMH=mk_meanMedian_bin(log10(hostmm),ncOut,'nbins',20);

gmass=galMass(galaxyMask);
ncGM=mk_meanMedian_bin(log10(gmass),ncGal,'nbins',20);
ncCM=mk_meanMedian_bin(log10(gmass),ncCGM,'nbins',20);
ncOM=mk_meanMedian_bin(log10(gmass),ncOut,'nbins',20);
%% all galaxies  

figure 

histogram(log10(ncGal),50);
hold on
histogram(log10(ncCGM),50);
histogram(log10(ncOut),50);

hl=legend('Gal','CGM','Out');
set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);

xlabelmine('log Number of non-SF gas cells');
ylabelmine('Number of galaxies');

titlemine(['All Galaxies ' simDisplayName])
set(gca,'Fontsize',14);
printout_fig(gcf,['cellNumHist_All_' simDisplayName]);
%% SF 

figure 

histogram(log10(ncGalSF),50);
hold on
histogram(log10(ncCGMSF),50);
histogram(log10(ncOutSF),50);

hl=legend('Gal','CGM','Out');
set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);

xlabelmine('log Number of non-SF gas cells');
ylabelmine('Number of galaxies');

titlemine(['Star-Forming Galaxies ' simDisplayName])
set(gca,'Fontsize',14);
printout_fig(gcf,['cellNumHist_SF_' simDisplayName]);
% Quenched 

figure 

histogram(log10(ncGalQ),50);
hold on
histogram(log10(ncCGMQ),50);
histogram(log10(ncOutQ),50);

hl=legend('Gal','CGM','Out');
set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);

xlabelmine('log Number of non-SF gas cells');
ylabelmine('Number of galaxies');

titlemine(['Quenched Galaxies ' simDisplayName])
set(gca,'Fontsize',14);
printout_fig(gcf,['cellNumHist_Q_' simDisplayName]);


%% host/galaxy mass dependence 

figure 

semilogy(ncGM.xMedian,ncGM.yMedian)
hold on
semilogy(ncCM.xMedian,ncCM.yMedian)
semilogy(ncOM.xMedian,ncOM.yMedian)
grid

hl=legend('Gal','CGM','Out');
set(hl,'Interpreter','latex','fontsize',16,'location','northwest');

xlabelmine('log Stellar Mass');
ylabelmine('mean cell number');

printout_fig(gcf,['cellNum_galMass_' simDisplayName]);

figure 

semilogy(ncGMH.xMedian,ncGMH.yMedian)
hold on
semilogy(ncCMH.xMedian,ncCMH.yMedian)
semilogy(ncOMH.xMedian,ncOMH.yMedian)
grid

hl=legend('Gal','CGM','Out');
set(hl,'Interpreter','latex','fontsize',16,'location','northwest');

xlabelmine('log Halo Mass');
ylabelmine('mean cell number');

printout_fig(gcf,['cellNum_haloMass_' simDisplayName]);


%% 
fname=sprintf('cellNumbers_snp%s_%s.mat',num2str(snap),simDisplayName);
    save([DEFAULT_MATFILE_DIR '/' fname],'cellNumStruct','-v7.3')
    
    fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
    



