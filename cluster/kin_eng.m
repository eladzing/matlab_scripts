function [Ekin Ekin_r] = kin_eng(boxx,vcmflag,rvcm,hubflag)
%creates cartesian cube of kinetic energy per volume
%returns values in solarmass * (km/sec)^2 per Mpc^3

%total kinetic energy
global hub

center=[0,0,0];

%load relevant values
rog=RHOG(boxx);

hf=1;
if exist('hubflag','var')
    hf=~strcmpi(hubflag,'nohub');    
end

if ~exist('rvcm','var')
    rvcm=[0 1];
end

if strcmp(vcmflag,'vcm') %%find Vrcm within Rvir
%     rvc=rvcm.*get_rvir();
    [VcmX VcmY VcmZ]=Vcm_full(rvcm,center,hubflag);
else
    global VCM
    VcmX=VCM(1); VcmY=VCM(2); VcmZ=VCM(3); % zeros(size(rog)); VcmY=zeros(size(rog)); VcmZ=zeros(size(rog));    
end
clear rvc;

[hubX hubY hubZ] = hubble_flow(boxx,center);
Vxx = Vx(boxx)+hf.*hubX-VcmX;
Vyy = Vy(boxx)+hf.*hubY-VcmY;
Vzz = Vz(boxx)+hf.*hubZ-VcmZ;

Vsq=Vxx.^2 +Vyy.^2 +Vzz.^2; %find v^2

clear VcmX VcmY VcmZ hubX hubY hubZ 

%% find Vr
lenc= size(rog);

[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));

%convert to center origin coordinates
meshX = meshX - (lenc(1)+1)/2 -center(1);
meshY = meshY - (lenc(2)+1)/2 -center(2);
meshZ = meshZ - (lenc(3)+1)/2 -center(3);
% Fix Units (to be in Mpc)

meshX = meshX * ((boxx/hub)/lenc(1));
meshY = meshY * ((boxx/hub)/lenc(2));
meshZ = meshZ * ((boxx/hub)/lenc(3));

rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2); % r^2 cube in Mpc

vr=(Vxx.*meshX+Vyy.*meshY+Vzz.*meshZ)./rcube;
clear Vxx Vyy Vzz 

Ekin=0.5.*rog.*Vsq;
Ekin_r=0.5.*rog.*vr.*vr;

% multiply by km / (1e6 * Mpc) * yr
%ff    = 0.5 .* Rhog.*(Vrr)^2.*ds.*(1.018120703e-12);