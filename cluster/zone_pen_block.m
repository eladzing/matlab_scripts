function zone_pen_block(plotproj,boxx,clims,hubbleflag)

%clustername='CL6';
%aexp='a1';typ='csf';
%new_env(clustername,typ,aexp);
%hubbleflag='hub';
%boxx=8;
%plotproj=[1 0 0];
zone=[0.01 0.2 1];
dilute = 8;
%clims=[-3 3];

%global CLUSTER;
global HALO_PATH;
global aexpn;
global NCELL;
global hub 

h=hub;
nthick=(1./NCELL).*4; 

lims=zeros([5 4 2]);

thick=nthick.*boxx;

load(sprintf('%s/virial%d_%s.mat', HALO_PATH, 8,aexpn));

[Mg200 ms200 M_dm200]=read_Mass_Profiles(RVIR);
    
%svir=interp1(R_Profile,S_Profile,RVIR,'spline');
%SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);

fluxnorm=(0.056.*(MVIR./1e13).^0.15)./(4.*pi);
ro=RHOG(boxx);
%tm=T(boxx);

%cson=sqrt(1.52e8.*tm)./1e5; %speed of sound in km/sec

%s=S(boxx)./SVIR;

%%find Vcm
rvfac=0.2;
rvcm=RVIR*rvfac;
vbox=ceil(2^ceil(log2(2*rvcm*h)));

rc=mk_rcube(vbox,ones(size(ro)));
rv_ind=find(rc<=rvcm);
[VcmX VcmY VcmZ] = V_Vcm_r(vbox,rv_ind,hubbleflag); 
clear rv_ind rc vbox;


if ~exist('cm')
    cm = [0,0,0];
end

hf=(strcmpi(hubbleflag,'hub'));    
[hubX hubY hubZ] = hubble_flow(boxx,cm);
vx = Vx(boxx)+hf.*hubX-VcmX;   
vy = Vy(boxx)+hf.*hubY-VcmY;   
vz = Vz(boxx)+hf.*hubZ-VcmZ;   
clear VcmX VcmY VcmZ hubX hubY hubZ;

%%Vtot=sqrt(Vx.^2+Vy.^2+Vz.^2);
[meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
%convert to center origin coordinates
meshX = meshX - (size(vx,1)+1)/2 -cm(1);
meshY = meshY - (size(vy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
%h=0.7;
meshX = meshX * ((boxx/h)/NCELL);
meshY = meshY * ((boxx/h)/NCELL);
meshZ = meshZ * ((boxx/h)/NCELL);

rcube2=meshX.^2+meshY.^2+meshZ.^2 ; % r^2 cube in Mpc
vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./sqrt(rcube2) ; %radial velocity 

hubflo=hub.*100.*sqrt(rcube2);
vrflo=vr+hf.*hubflo;

flux=ro.*(vrflo).*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ; 
%flux=ro.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ; 

clear meshY meshX meshZ rcube2 vrflo hubflo;


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
%resampled_v_x = imresize(v_x, [diluted_len diluted_len], 'bilinear');
%resampled_v_y = imresize(v_y, [diluted_len diluted_len], 'bilinear');

tickjump = boxx/4;
load('MyColormaps','avijet');

%% plot flux

%texttowrite='$dot M /M_{\mathrm{vir}} [1/\mathrm{Gyr}]$';
%clims=[-3 3];
mappin_flux

set(gca,...
    'Position',[0.08485 0.124 0.7844 0.8557],...
    'FontSize',14,...
    'DataAspectRatio',[1 1 1]);%'PlotBoxAspectRatio',[1 1 1],...
   % 'Position',[0.07 0.55 0.37 0.37],...
    bar=colorbar('position',[0.8419 0.1196 0.05051 0.8193]);
    %'Position',[0.1634 0.1316 0.7079 0.8519],...
    %[0.8929 0.1345 0.06188 0.8363]);
 set(get(bar,'Title'),'String','$\frac{\dot M}{d\Omega} [\frac{\dot M}{d\Omega}\big |_{\mathrm{virial}}]$','Fontsize',14,'Interpreter','latex');

 
 
% ylabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');
% xlabel('$\mathrm{[Mpc/h]}$','Fontsize',12,'Interpreter','latex');

