%% plottting script for 

%% load data 
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/cooling_times_z0_TNG100.mat')
load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/tng100_z0_fofs_subs.mat')

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
units; % load general unit structure in cgs.

% identiy central s 
centralMask= subsInfo.isCentral(tCoolStruct.galMask);

% get ssfr 
gmass=tCoolStruct.galMass(tCoolStruct.galMask);
ssfr=subs.SubhaloSFRinRad(tCoolStruct.galMask)./gmass+1e-16;
tc=tCoolStruct.inCGM.meanTcMW(:,1)';

%% plot cooling time vs. mass gal overlaid on mass -ssfr for all
% filt=fspecial('disk',4);
% res=plot_population_contour(log10(gmass),log10(tc'),...
% 'xlim',[9 13],'ylim',[-6 2],'smooth',filt);

xlabelmine('$\log M_\star\,[\mathrm{M_\odot}]$')
ylabelmine('$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$')
%titlemine(['All galaxies:' num2str(length(mass))])

%% plot cooling time vs. mass gal overlaid on mass -ssfr for all
figure

cmap=brewermap(256,'*Spectral');
scatter(log10(gmass(centralMask)),log10(ssfr(centralMask)),10,log10(tc'))
colormap(cmap)
hb=colorbar;barTitle(hb,'$\log t_{\mathrm{cool}}\,[\mathrm{Gyr}]$')
xlabelmine('$\log M_\star\,[\mathrm{M_\odot}]$')
ylabelmine('$\log \mathrm{ssfr}\,[\mathrm{yr^{-1}}] $') ;
%% plot cooling time vs. mass gal overlaid on mass -ssfr for all
figure

%cmap=brewermap(256,'*Spectral');
scatter(log10(gmass(centralMask)),log10(tc(centralMask)),10,log10(ssfr(centralMask)))
colormap(cmap)
hb=colorbar;
xlabelmine('$\log M_\star\,[\mathrm{M_\odot}]$')
ylabelmine('$\log t_\mathrm{cool}\,[\mathrm{Gyr}]$')
barTitle(hb,'$\log \mathrm{ssfr}\,[\mathrm{yr^{-1}}] $')