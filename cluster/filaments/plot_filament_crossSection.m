function plot_filament_crossSection( csSt,varargin )
%PLOT_FILAMENT_CROSSSECTION - plot the filament cross section  
%   Functuib accepts the data strucutre containging the 
%   cross-section data. 

plotFlag=false(12,1);
if isempty(varargin)
    plotFlag(:)=true;
else
    i=1;
    while i<=length(varargin)
        switch(lower(varargin{i}))
            case 'all'
                plotFlag(:)=true;
            case {'bird','phase','rot'}
                plotFlag(1)=true;
            case {'ro','density'}
                plotFlag(2)=true;
            case {'t','tmp','temperature'}
                plotFlag(3)=true;
            case {'ent','entropy','s','k'}
                plotFlag(4)=true;
            case {'p','press','pre','pressure'}
                plotFlag(5)=true;
            case {'m','mach'}
                plotFlag(6)=true;
            case {'vt'}
                plotFlag(7)=true;
            case {'vr'}
                plotFlag(8)=true;
            case {'mu'}
                plotFlag(9)=true;
            case {'flux','f'}
                plotFlag(10)=true;
            case {'vp'}
                plotFlag(11)=true;
            case {'dm'}
                plotFlag(12)=true;
            otherwise
                error('plot_filament_crossSection - Illegal argument: %s',varargin{i});
        end
        i=i+1;
    end
end

%% plot
pTag=csSt.pTag;
printoutdir=csSt.printoutDir;
dPl=csSt.dPLen;
xLab=sprintf('$t_p\\,[%s\\,\\mathrm{kpc}]$',num2str(dPl));


% bird
if plotFlag(1)
figure
scatter(log10(csSt.ro),log10(csSt.tmp),15,csSt.tParam)
cmap=brewermap(256,'*Spectral');
colormap(cmap)
set(gca,'color',[0.6 0.6 0.6])
bar=colorbar;
set(bar,'fontsize',10);
barTitle(bar,xLab)
caxis([0 1])
xlabelmine('$\log n_{gas}\,[\mathrm{cm^{-3}}]$');
ylabelmine('$\log T\,[\mathrm{K}]$');
grid
set(gca,'box','on','Fontsize',14);
fname=sprintf('%s_filCS_phase_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end



% density
if plotFlag(2)
figure
semilogy(csSt.tParam,csSt.ro,'.')
xlabelmine(xLab)
ylabelmine('$\log n_{gas}\,[\mathrm{cm^{-3}}]$');
grid
fname=sprintf('%s_filCS_rho_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% temperature
if plotFlag(3)
figure
semilogy(csSt.tParam,csSt.tmp,'.')
xlabelmine(xLab)
ylabelmine('$\log T\,[\mathrm{K}]$');
grid
fname=sprintf('%s_filCS_tmp_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% entropy
if plotFlag(4)
figure
semilogy(csSt.tParam,csSt.ent,'.')
xlabelmine(xLab)
ylabelmine('$\log S\,[\mathrm{KeV\,cm^2}]$');
grid
fname=sprintf('%s_filCS_ent_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% pressure
if plotFlag(5)
figure
semilogy(csSt.tParam,csSt.pres,'.')
xlabelmine(xLab)
ylabelmine('$\log P\,[\mathrm{dym\,cm^{-2}}]$');
grid
fname=sprintf('%s_filCS_press_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% mach
if plotFlag(6)
figure
semilogy(csSt.tParam,csSt.mach,'.')
xlabelmine(xLab)
ylabelmine('local $\mathcal{M}$');
grid
fname=sprintf('%s_filCS_mach_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% vtot
if plotFlag(7)
figure
plot(csSt.tParam,csSt.vt,'.')
xlabelmine(xLab)
ylabelmine('$v/V_{\mathrm{vir}}$');
grid
fname=sprintf('%s_filCS_vtot_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% vr
if plotFlag(8)
figure
plot(csSt.tParam,csSt.vr,'.')
xlabelmine(xLab)
ylabelmine('$v_r/V_{\mathrm{vir}}$');
grid
fname=sprintf('%s_filCS_vr_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% mu
if plotFlag(9)
figure
plot(csSt.tParam,csSt.mu,'.')
xlabelmine(xLab)
ylabelmine('$\mu$');
grid
fname=sprintf('%s_filCS_mu_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% flux
if plotFlag(10)
figure
plot(csSt.tParam,csSt.flux,'.')
xlabelmine(xLab)
ylabelmine('$\dot{M}_{gas}/\mathrm{d}\Omega$ [arbit. units]');
grid
fname=sprintf('%s_filCS_flux_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% perpindicular velocity
if plotFlag(11)
figure
plot(csSt.tParam,csSt.vp,'.')
xlabelmine(xLab)
ylabelmine('$v_{\perp}/v_{\mathrm{vir}}$');
grid
fname=sprintf('%s_filCS_vperp_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end

% dm
if plotFlag(12)
figure
semilogy(csSt.tParam,csSt.dm,'.')
xlabelmine(xLab)
ylabelmine('$\rho_{\mathrm{DM}}\,[\mathrm{M_\odot \, Mpc^{-3}}]$')
grid
fname=sprintf('%s_filCS_dm_cs%s',csSt.cluster,num2str(pTag));
printout_fig(gcf,fname,'dir',printoutdir,'png');
end


end

