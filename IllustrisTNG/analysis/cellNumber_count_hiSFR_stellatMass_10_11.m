%% plot the number of gas cells (without star-forming gas) 
%  within the Gal, CGM, Out, components, also seperating into SF and Quenched galaxies 
snap=99;
global simDisplayName
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

massMask=galMass>1e10 & galMass<1e11;

ssfrThresh=10^(-9.5);
qMask=ssfr<=ssfrThresh;

galaxyMask=centralMask & ~spbMask(tCoolStruct.galMask) & massMask;

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
%printout_fig(gcf,['cellNumHist_All_' simDisplayName]);
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

titlemine(['Star-Forming Galaxies above -9.5, transition mass ' simDisplayName])
set(gca,'Fontsize',14);
printout_fig(gcf,['cellNumHist_hiSF_' simDisplayName]);
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
%printout_fig(gcf,['cellNumHist_Q_' simDisplayName]);

% 
% %% 
% fname=sprintf('cellNumbers_snp%s_%s.mat',num2str(snap),simDisplayName);
%     save([DEFAULT_MATFILE_DIR '/' fname],'cellNumStruct','-v7.3')
%     
%     fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);
%     



