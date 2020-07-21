zone=[0.01 0.2 1];
nthick=(1./256).*4; 

%lims=zeros([5 4 2]);
%clustername='CL6';
new_env(clustername);
%plotproj=[1 0 0];

global FILE_FORMAT;
global HALO_PATH
load(sprintf('%s/virial%d', HALO_PATH, 8));

    
svir=interp1(R_Profile,S_Profile,RVIR,'spline');
SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);

fluxnorm=(0.056.*(MVIR./1e13).^0.15)./(4.*pi);

figure;
%% plot 4 box 

boxx=4;

thick=nthick.*boxx;
cson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec
ro=RHOG(boxx);
tm=T(boxx);
s=S(boxx)./SVIR;

%%find Vcm
h=0.7;
rvfac=0.2;
rvcm=RVIR*rvfac;
vbox=ceil(2^ceil(log2(2*rvcm*h)));

% switch clustername
%     case {'CL101','CL102','CL103','CL104','CL105','CL106','CL107','CL5'}
%         if vbox==1
%             vbox=2;
%         end
% end

rc=mk_rcube(vbox,ones(size(ro)));
rv_ind=find(rc<=rvcm);
[VcmX VcmY VcmZ] = V_Vcm_r(vbox,rv_ind); 
clear rv_ind rc vbox;

if ~exist('cm')
    cm = [0,0,0];
end

vx = Vx(boxx)-VcmX;   vy = Vy(boxx)-VcmY;   vz = Vz(boxx)-VcmZ;   
clear VcmX VcmY VcmZ;

%%Vtot=sqrt(Vx.^2+Vy.^2+Vz.^2);
[meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
%convert to center origin coordinates
meshX = meshX - (size(vx,1)+1)/2 -cm(1);
meshY = meshY - (size(vy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
%h=0.7;
meshX = meshX * ((boxx/h)/256);
meshY = meshY * ((boxx/h)/256);
meshZ = meshZ * ((boxx/h)/256);

rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity 

flux=ro.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ; 

Mach=abs(vr./cson); clear vr cson;
rho=ro./(3.*MVIR./(4.*pi.*RVIR.^3));

clear meshY meshX meshZ rcube2;

%% plot maps

dilute = 8;
blen=size(ro,1);
thk=ceil(0.5.*thick./boxx.*blen);   %% thick is in Mpc
slind=(floor(0.5.*blen)-thk+1):1:(ceil(0.5.*blen)+thk);
side=[-boxx./2 boxx./2];

zones=zone(find(zone<=(0.5.*boxx./h./RVIR)));

[v_x v_y]=mk_vfield(vx,vy,vz,ro,slind); 
clear vx vy vz;
diluted_len = length(1:dilute:blen);
diluted_jump = boxx/(diluted_len-1);
notdiluted_jump = boxx/(blen-1);
[xxv yyv] = meshgrid(-boxx/2:diluted_jump:boxx/2, -boxx/2:diluted_jump:boxx/2);
[xxs yys] = meshgrid(-boxx/2:notdiluted_jump:boxx/2,-boxx/2:notdiluted_jump:boxx/2);
resampled_v_x = imresize(v_x, [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(v_y, [diluted_len diluted_len], 'bilinear');
tickjump = boxx/4;
load('MyColormaps','avijet');

%% plot flux
subplot(2,2,1);
%texttowrite='$dot M /M_{\mathrm{vir}} [1/\mathrm{Gyr}]$';
clims=[-3 3];
mappin_flux

set(gca,...
    'Position',[0.07 0.55 0.37 0.37],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',12,...
    'DataAspectRatio',[1 1 1]);
 bar=colorbar('position',[0.438 0.55 0.033 0.37]);
 set(get(bar,'Title'),'String','$\frac{\dot M/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\big |_{\mathrm{virial}}]$','Fontsize',12,'Interpreter','latex');
ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');

%% plot Mach
subplot(2,2,2);
clims=[-2 2];
%texttowrite='$T [\mathrm{K}]$';
mappin_mach; %    zoomapt;
 set(gca,...
    'Position',[0.55 0.55 0.37 0.37],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',12,...
    'DataAspectRatio',[1 1 1]);
 bar=colorbar('position',[0.918 0.55 0.033 0.37]);
 set(get(bar,'Title'),'String','$\mathrm{log}\left(|v_r/c_{s}|\right)$','Fontsize',12,'Interpreter','latex');

%ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
%xlabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');

%% plot 1 box 
clear boxx thick rho ro tm s Mach flux v_x v_y diluted_len diluted_jump notdiluted_jump;
clear h xxv yyv xxs yys resampled_v_x resampled_v_y slice;

boxx=1;

thick=nthick.*boxx;
cson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec
ro=RHOG(boxx);
tm=T(boxx);
s=S(boxx)./SVIR;

%%find Vcm
h=0.7;
rvfac=0.2;
rvcm=RVIR*rvfac;
vbox=ceil(2^ceil(log2(2*rvcm*h)));

switch clustername
    case {'CL101','CL102','CL103','CL104','CL105','CL106','CL107','CL5'}
        if vbox==1
            vbox=2;
        end
end

rc=mk_rcube(vbox,ones(size(ro)));
rv_ind=find(rc<=rvcm);
[VcmX VcmY VcmZ] = V_Vcm_r(vbox,rv_ind); 
clear rv_ind rc vbox;

if ~exist('cm')
    cm = [0,0,0];
end

vx = Vx(boxx)-VcmX;   vy = Vy(boxx)-VcmY;   vz = Vz(boxx)-VcmZ;   
clear VcmX VcmY VcmZ;

%%Vtot=sqrt(Vx.^2+Vy.^2+Vz.^2);
[meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
%convert to center origin coordinates
meshX = meshX - (size(vx,1)+1)/2 -cm(1);
meshY = meshY - (size(vy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
%h=0.7;
meshX = meshX * ((boxx/h)/256);
meshY = meshY * ((boxx/h)/256);
meshZ = meshZ * ((boxx/h)/256);

rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity 

flux=ro.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ; 

Mach=abs(vr./cson); clear vr cson;
rho=ro./(3.*MVIR./(4.*pi.*RVIR.^3));

clear meshY meshX meshZ rcube2;

%% plot maps

dilute = 8;
blen=size(ro,1);
thk=ceil(0.5.*thick./boxx.*blen);   %% thick is in Mpc
slind=(floor(0.5.*blen)-thk+1):1:(ceil(0.5.*blen)+thk);
side=[-boxx./2 boxx./2];

zones=zone(find(zone<=(0.5.*boxx./h./RVIR)));

[v_x v_y]=mk_vfield(vx,vy,vz,ro,slind); 
clear vx vy vz;
diluted_len = length(1:dilute:blen);
diluted_jump = boxx/(diluted_len-1);
notdiluted_jump = boxx/(blen-1);
[xxv yyv] = meshgrid(-boxx/2:diluted_jump:boxx/2, -boxx/2:diluted_jump:boxx/2);
[xxs yys] = meshgrid(-boxx/2:notdiluted_jump:boxx/2,-boxx/2:notdiluted_jump:boxx/2);
resampled_v_x = imresize(v_x, [diluted_len diluted_len], 'bilinear');
resampled_v_y = imresize(v_y, [diluted_len diluted_len], 'bilinear');
tickjump = boxx/4;
load('MyColormaps','avijet');


%% plot flux
   subplot(2,2,3);
   clims=[-3 3];
   %texttowrite='$\rho/\rho_{\mathrm{vir}}$';
   mappin_flux;
   set(gca,...
    'Position',[0.07 0.09 0.37 0.37],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',12,...
    'DataAspectRatio',[1 1 1]);
   bar=colorbar('position',[0.438 0.09 0.033 0.37]);
  %set(get(bar,'Title'),'String','$\frac{\dot M/M_{gas}}{d\Omega} [\frac{\dot M}{M_{gas}}\big |_{\mathrm{virial}}]$','Fontsize',12,'Interpreter','latex');
ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
xlabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');

  
%% plot Mach
subplot(2,2,4);
%texttowrite='$S/S(R_{\mathrm{vir}})$';
clims=[-2 2];
mappin_mach;
 set(gca,...
    'Position',[0.55 0.09 0.37 0.37],...
    'PlotBoxAspectRatio',[1 1 1],...
    'FontSize',12,...
    'DataAspectRatio',[1 1 1]);
    
bar=colorbar('position',[0.918 0.09 0.033 0.37]);
%set(get(bar,'Title'),'String','$\mathrm{log}\left(\big|\frac{v_r}{c_{s}}\big|\right)$','Fontsize',12,'Interpreter','latex');
xlabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
