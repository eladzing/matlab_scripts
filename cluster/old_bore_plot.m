%% read 
%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
cl=6;cltype='csf';aexp='a1';
pflag='noprint';

clustername=sprintf('CL%d',cl)
new_env(clustername,cltype,aexp);

boxx=1;

global HALO_PATH
load(sprintf('%s/virial%d', HALO_PATH, 8));
    
%vir=interp1(R_Profile,S_Profile,RVIR,'spline');





SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3);

fluxnorm=(0.056.*(MVIR./1e13).^0.15)./(4.*pi);

cson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec


ro=RHOG(boxx);
tm=T(boxx);
s=S(boxx)./SVIR;
cson=sqrt(1.52e8.*T(boxx))./1e5; %speed of sound in km/sec
pre=(ro.*tm)/(TVIR.*(3.*MVIR./(4.*pi.*RVIR.^3)));

%%find Vcm
h=0.7;
rvfac=0.2;
rvcm=RVIR*rvfac;
vbox=ceil(2^ceil(log2(2*rvcm*h)));

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

%flux=ro.*vr.*rcube2./Mg200.*(1e5./3.0856e24.*3.15e7.*1e9)./fluxnorm ; 

%Mach=abs(vr./cson); clear vr cson;

[gsx gsy gsz]=gradient(s);

grads=sqrt(gsx.^2+gsy.^2+gsz.^2);   
gradsrad=(gsx.*meshX+gsy.*meshY+gsz.*meshZ)./sqrt(rcube2) ;

clear gsx gsy gsz meshY meshX meshZ rcube2;


%% plot maps

%dilute = 4;
blen=size(ro,1);
thk=ceil(0.5.*thick./boxx.*blen);   %% thick is in Mpc
slind=(floor(0.5.*blen)-thk+1):1:(ceil(0.5.*blen)+thk);

if cshift~=0
    shift_cen=0.5.*blen-(ciel(0.5.*blen*(cshift/boxx*2+1)));
else
    shift_cen=0;
end

if zmbox~=0
    zmind=ceil(zmbox./boxx.*blen);
else
    zmind=0.5.*blen;
end

slind=slind-shift_cen;
sidind=(0.5*blen+1)-zmind:0.5*blen+zmind;
actbox=zmind*(boxx/blen);
side=[-actbox actbox];
zblen=length(sidind);

zones=zone(find(zone<=(0.5.*boxx./h./RVIR)));


[v_x v_y]=mk_vfield(vx,vy,vz,ro,slind); 
clear vx vy vz;
diluted_len = length(1:dilute:zblen);
diluted_jump = 2.*actbox/(diluted_len-1);
notdiluted_jump = 2.*actbox/(zblen-1);
[xxv yyv] = meshgrid(-actbox:diluted_jump:actbox, -actbox:diluted_jump:actbox);
[xxs yys] = meshgrid(-actbox:notdiluted_jump:actbox,-actbox:notdiluted_jump:actbox);
%resampled_v_x = imresize(v_x, [diluted_len diluted_len], 'bilinear');
%resampled_v_y = imresize(v_y, [diluted_len diluted_len], 'bilinear');

tickjump = actbox/3;
load('MyColormaps','avijet');


%% plot entropy
% printag='s_lin';
% clims=[-1.5 0.5]; %  s_lim';
% texttowrite='(S/S(r=R_{vir}))';
% if strcmp(tag3,'yes')
%     mappin3_ent;
% else
%     mappin_ent;
%     %clims=[-2 0.25];
%     %zoomaps;
% end

%% plot temp
% printag='tm';
% clims=[6 8.5];%t_lim';
% %texttowrite='(T/T_vir)';
% if strcmp(tag3,'yes')
%     mappin3_tm;
% else
%     mappin_tm;
% end

%% plot pressure
printag='pre';
clims=[-2 3];%t_lim';
texttowrite='(T/T_vir)';
if strcmp(tag3,'yes')
    mappin3_pre;
else
    mappin_pre;      
end


% %% plot mach
% 
% printag='mach_vzone2';
% clims=m_lim';
% texttowrite='log(|v_r/c_s|)';
% 
% % if strcmp(tag3,'yes')
% %     mappin3_mach;
% % else
% %     mappin_mach
% % end
% 
% %% plot flux
% printag='flux_vzone2';
% clims=f_lim';
% texttowrite='Mdot/Mdot_{vir}]';
%  if strcmp(tag3,'yes')
%      mappin3_flux;
%  else
%      clims=[-3 3];
%      zoomapfl;
%      %mappin_flux
%  end
% 
%% plot entropy gradient
% [gsx gsy gsz]=gradient(s);
% 
% grads=sqrt(gsx.^2+gsy.^2+gsz.^2);clear gsx gsy gsz;
% 
% clims=[0 0]; %gs_lim';
% printag='grad_s';
% texttowrite='log(\nabla S/S_{vir})';
% 
% if strcmp(tag3,'yes')
%     mappin3_gradent;
% else
%     mappin_gradent 
% end

% 
% %% plot radial component of entropy gradient
% 
% % clims=gsr_lim';
% % printag='grad_s_rad';
% % texttowrite='log(\nabla S/S_{vir})';
% % 
% % if strcmp(tag3,'yes')
% %     %mappin3_gradentrad;
% % else
% %     %mappin_gradentrad; 
% % end
% 
