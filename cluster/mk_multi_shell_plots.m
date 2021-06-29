%% plot hammer plots for mass inflow, density temperature and entropy
cl=9;
basename='Cl%s_%s_hammer_r%s';
units;
new_env(cl,'a06');

global NCELL
global hub
global zred


thik= 2 ;%thickness of shell

rotThet=3*pi/2;
rotPhi=0; % should always be zero since the sampleing is dependent on this!
mv=get_mvir;
rv=get_rvir;
tv=get_tvir;

map=brewermap(256,'*RdYlBu');
%load('MyColormaps','newJet');

rom=mv/(4*pi/3*rv^3);
fluxnorm=(0.056.*(mv./1e13).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
%% identify radius
rad=0.26;%   0.7/hub/rv;  %in units of Rvir

radh=rad*rv.*hub;

boxx=2.^(max(ceil(log(radh)/log(2))+1,0));

if ~any(boxx==[1 2 4 8])
    error('mk_multi_shell_plots: Illegal box %i',boxx);
end
[~,mesh_phi,mesh_theta] = sphere_grid(boxx);


indx=ceil(radh/(0.5*boxx)*NCELL);

ll=indx-thik:indx+thik;

wt=RHOG_sphere(boxx);

%% mass inflow
cu=flux_sphere(boxx)./fluxnorm;
shell=squeeze(sum(cu(ll,:,:),1))/length(ll);

mphi=squeeze(mesh_phi(indx,:,:));
mthet=squeeze(mesh_theta(indx,:,:));
%rad=ind/NCELL*0.5*boxx/(1+zred)/hub/get_rvir;

% rotate

dThet=abs(max(max(diff(mthet,1,2))));
%dPhi=abs(max(max(diff(mphi,1,1))));

if rotThet~=0
    dInd=floor(rotThet/dThet);
    sh=zeros(size(shell));
    sh(:,dInd+1:end)=shell(:,1:end-dInd);
    sh(:,1:dInd)=shell(:,end-dInd+1:end);
    shell=sh;
    clear sh;
end
% rotation along this variable left only for educational purposes  
% if rotPhi~=0
%     dInd=floor(rotPhi/dPhi);
%     sh=zeros(size(shell));
%     sh(dInd+1:end,:)=shell(1:end-dInd,:);
%     sh(1:dInd,:)=shell(end-dInd+1:end,:);
%     shell=sh;
%     clear sh;
% end
% 

hf=figure;
plot_hammer_shell(shell,mthet,mphi,'nc',500);
caxis([-30 30])
colormap(map);  
bar=colorbar;set(bar,'fontsize',14)
barTitle(bar,'$\frac{\dot{M}}{d\Omega}\,[\frac{\dot{M}}{d\Omega}\big|_{\mathrm{vir}}]$');
%titlemine(sprintf('CL%s',num2str(cl)));
name=sprintf(basename,num2str(cl),'flux',num2str(rad));
set(gcf,'Position',[402   566   718   382]);
%saveas(hf,sprintf('%s.fig',name));
% printout _fig(hf,name,'v')

%% entropy

cu=S_sphere(boxx).*f_ent; %/(kb*tv/1e3/ev);
shell=squeeze(sum(cu(ll,:,:).*wt(ll,:,:),1)./sum(wt(ll,:,:),1));

mphi=squeeze(mesh_phi(indx,:,:));
mthet=squeeze(mesh_theta(indx,:,:));
%rad=ind/NCELL*0.5*boxx/(1+zred)/hub/get_rvir;

%SVIR=tv.*rv.^2./(3.*get_mvir./(4.*pi)).^(2./3);

% rotate

dThet=abs(max(max(diff(mthet,1,2))));
%dPhi=abs(max(max(diff(mphi,1,1))));

if rotThet~=0
    dInd=floor(rotThet/dThet);
    sh=zeros(size(shell));
    sh(:,dInd+1:end)=shell(:,1:end-dInd);
    sh(:,1:dInd)=shell(:,end-dInd+1:end);
    shell=sh;
    clear sh;
end

hf=figure;
plot_hammer_shell(log10(shell),mthet,mphi,'nc',500);
caxis([2 3.1])
colormap(map);  
bar=colorbar;set(bar,'fontsize',14)
barTitle(bar,'$\log S \,[\mathrm{KeV\,cm^2}]$');
%titlemine(sprintf('CL%s',num2str(cl)));
name=sprintf(basename,num2str(cl),'entropy',num2str(rad));
%saveas(hf,sprintf('%s.fig',name));
set(gcf,'Position',[402   566   718   382]);
% printout _fig(hf,name,'v')
%% density
cu=RHOG_sphere(boxx);
shell=squeeze(sum(cu(ll,:,:),1))./length(ll)./rom;

mphi=squeeze(mesh_phi(indx,:,:));
mthet=squeeze(mesh_theta(indx,:,:));
%rad=ind/NCELL*0.5*boxx/(1+zred)/hub/get_rvir;

% rotate

dThet=abs(max(max(diff(mthet,1,2))));
%dPhi=abs(max(max(diff(mphi,1,1))));

if rotThet~=0
    dInd=floor(rotThet/dThet);
    sh=zeros(size(shell));
    sh(:,dInd+1:end)=shell(:,1:end-dInd);
    sh(:,1:dInd)=shell(:,end-dInd+1:end);
    shell=sh;
    clear sh;
end

hf=figure;
plot_hammer_shell(shell,mthet,mphi,'nc',500);
caxis([0.2 0.7])
colormap(map);  
bar=colorbar;
set(bar,'fontsize',14)
barTitle(bar,'$\rho_{\mathrm{gas}}/\rho_{\mathrm{vir}}$');
%titlemine(sprintf('CL%s',num2str(cl)));
name=sprintf(basename,num2str(cl),'density',num2str(rad));
%saveas(hf,sprintf('%s.fig',name));
set(gcf,'Position',[402   566   718   382]);
% printout _fig(hf,name,'v')

%% temperature
cu=T_sphere(boxx)./tv;
shell=squeeze(sum(cu(ll,:,:).*wt(ll,:,:),1)./sum(wt(ll,:,:),1));
%shell=squeeze(sum(cu(ll,:,:),1))/length(ll);

mphi=squeeze(mesh_phi(indx,:,:));
mthet=squeeze(mesh_theta(indx,:,:));
%rad=ind/NCELL*0.5*boxx/(1+zred)/hub/get_rvir;

% rotate

dThet=abs(max(max(diff(mthet,1,2))));
%dPhi=abs(max(max(diff(mphi,1,1))));

if rotThet~=0
    dInd=floor(rotThet/dThet);
    sh=zeros(size(shell));
    sh(:,dInd+1:end)=shell(:,1:end-dInd);
    sh(:,1:dInd)=shell(:,end-dInd+1:end);
    shell=sh;
    clear sh;
end

hf=figure;
plot_hammer_shell(shell,mthet,mphi,'nc',500);
caxis([0.2 1.7])
colormap(map);  
bar=colorbar;set(bar,'fontsize',14)
barTitle(bar,'$T/T_{\mathrm{vir}}$');
%titlemine(sprintf('CL%s',num2str(cl)));
name=sprintf(basename,num2str(cl),'temperature',num2str(rad));
%saveas(hf,sprintf('%s.fig',name));
set(gcf,'Position',[402   566   718   382]);
% printout _fig(hf,name,'v')