%% Plot results for gas properties in Centrals in TNG
%% load data
%%
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/cooling_times_z0_TNG100.mat')
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/tng100_z0_fofs_subs.mat')

fofs=illustris.utils.addTvirFofs(fofs);


subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.

% identify centrals 
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

% get usefel stuff 
galMass=tCoolStruct.galMass(tCoolStruct.galMask)';  % galaxy stellar mass
gasMass=tCoolStruct.inGal.mass;
sfr=subs.SubhaloSFRinRad(tCoolStruct.galMask)';  % sfr in galaxy 
ssfr=sfr./galMass+1e-16;
sfre(:,1)=sfr./gasMass(:,1);
sfre(:,2)=sfr./gasMass(:,2);



tc=tCoolStruct.inCGM.meanTcMW(:,1)';
gasTemp=tCoolStruct.inCGM.meanTempMW(:,1);
gasEnt=tCoolStruct.inCGM.meanEntMW(:,1);

% get tvir for host fofs 
tvir=fofs.Group_T_Mean200(subsInfo.hostFof+1);
tvir=tvir(tCoolStruct.galMask);
tvirMean=mk_meanMedian_bin(log10(galMass(centralMask)),log10(tvir(centralMask)),'nb',20);



%% plotting stuff
mstarLab='$\log M_\star\,[\mathrm{M_\odot}]$';
tcLab='$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$';
ssfrLab='$\log \mathrm{sSFR}\,[\mathrm{yr^{-1}}] $';
ssfeLab='$\log \mathrm{SFRe}\,[\mathrm{yr^{-1}}] $';
mgasLab='$\log M_\mathrm{gas}\,[\mathrm{M_\odot}]$';
tempLab='$\log T \,[\mathrm{K}]$';
entLab='$\log S \,[\mathrm{KeV\,cm^2}]$';
fgsLab='$\log M_\mathrm{gas}/M_\mathrm{star}$';
%cmap=brewermap(256,'RdYlBu');
cmap=brewermap(256,'Spectral');

%% Mass vs Temp


% ssfr
%figure
% subplot(2,2,1)
% scatter(log10(galMass(centralMask)),log10(gasTemp(centralMask)),10,log10(ssfr(centralMask)));
% hold on 
% errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(tempLab);
% hb=colorbar;barTitle(hb,ssfrLab)
% colormap(cmap)
% %printout_fig(gcf,
% 

% sfre
%figure
% subplot(2,2,2)
% scatter(log10(galMass(centralMask)),log10(gasTemp(centralMask)),10,log10(sfre(centralMask,2)));
% hold on 
% errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(tempLab);
% hb=colorbar;barTitle(hb,ssfeLab)
% colormap(cmap)

% tcool
figure
%subplot(2,2,3)
scatter(log10(galMass(centralMask)),log10(gasTemp(centralMask)),5,log10(tc(centralMask)));
hold on 
errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
grid
xlabelmine(mstarLab);
ylabelmine(tempLab);
hb=colorbar;barTitle(hb,tcLab)
colormap(cmap)
colormap(flipud(cmap))



% gasMass
%figure
% subplot(2,2,4)
% fgs=tCoolStruct.inGal.mass(:,2)./galMass;
% scatter(log10(galMass(centralMask)),log10(gasTemp(centralMask)),5,log10(fgs(centralMask)));
% hold on 
% errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(tempLab);
% hb=colorbar;barTitle(hb,fgsLab)
% %colormap(flipud(cmap))
% colormap(cmap)


% entropy
% figure
% scatter(log10(galMass(centralMask)),log10(gasTemp(centralMask)),5,log10(gasEnt(centralMask)));
% hold on 
% errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(tempLab);
% hb=colorbar;barTitle(hb,entLab)
% colormap(flipud(cmap))


%% mass vs entropy 
% % ssfr
% figure
% subplot(2,2,1)
% scatter(log10(galMass(centralMask)),log10(gasEnt(centralMask)),10,log10(ssfr(centralMask)));
% % hold on 
% % errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(entLab)
% hb=colorbar;barTitle(hb,ssfrLab)
% colormap(cmap)
% %printout_fig(gcf,
% 
% 
% % sfre
% %figure
% subplot(2,2,2)
% scatter(log10(galMass(centralMask)),log10(gasEnt(centralMask)),10,log10(sfre(centralMask,2)));
% % hold on 
% % errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(entLab)
% hb=colorbar;barTitle(hb,ssfeLab)
% colormap(cmap)
% 

% tcool
figure
%subplot(2,2,3)
scatter(log10(galMass(centralMask)),log10(gasEnt(centralMask)),5,log10(tc(centralMask)));
% hold on 
% errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
grid
xlabelmine(mstarLab);
ylabelmine(entLab)
hb=colorbar;barTitle(hb,tcLab)
%colormap(cmap)
colormap(flipud(cmap))



% % gasMass
% %figure
% subplot(2,2,4)
% fgs=tCoolStruct.inGal.mass(:,2)./galMass;
% scatter(log10(galMass(centralMask)),log10(gasEnt(centralMask)),5,log10(fgs(centralMask)));
% % hold on 
% % errorbar(tvirMean.bins,tvirMean.Mean,tvirMean.stanDev/2,'dk')
% grid
% xlabelmine(mstarLab);
% ylabelmine(entLab)
% hb=colorbar;barTitle(hb,fgsLab)
% %colormap(flipud(cmap))
% colormap(cmap)












