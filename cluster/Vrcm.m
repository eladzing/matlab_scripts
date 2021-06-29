

function vrcm=Vrcm(boxx,rvcm,cm,hubflag) 
%%returns cube of radial component of Vcm. Vcm is  calculated
%within a given radial zone.
if ~exist('cm')
    cm = [0,0,0];
end

if ~exist('hubflag')
    hf=0;
else
    hf=(strcmpi(hubflag,'hub'));    
end

%[hubX hubY hubZ] = hubble_flow(boxx,cm); -> probably an error


%% calculate Vcm within a given radius

%read velocity and density from appropriate box 
if length(rvcm)<2
    rvcm=[0 rvcm];
end

global hub
h=hub;
vbox=ceil(2^ceil(log2(2*max(rvcm(:))*h)));

[hubX hubY hubZ] = hubble_flow(vbox,cm);

Vxx = Vx(vbox)+hf.*hubX;
Vyy = Vy(vbox)+hf.*hubY;
Vzz = Vz(vbox)+hf.*hubZ;

Rhog = RHOG(vbox);
lenc= size(Rhog);

%find index list of the volume 
rc=mk_rcube(vbox,ones(size(Rhog)),cm);
rind=find(rc>=min(rvcm(:)) & rc<=max(rvcm(:)));

%calculate Vcm
Rhog = Rhog(rind);
M = sum(Rhog(:));

vxro=Vxx(rind).*Rhog;
vyro=Vyy(rind).*Rhog;
vzro=Vzz(rind).*Rhog;

VcmX = sum(vxro(:))/M;
VcmY = sum(vyro(:))/M;
VcmZ = sum(vzro(:))/M;



%% prepare cube.
[meshY, meshX, meshZ] = meshgrid(1:size(Vyy,1), 1:size(Vxx,2), 1:size(Vzz,3));
clear Vxx Vyy Vzz Rhog vxro vyro vzro rind;


%convert to center origin coordinates
meshX = meshX - (lenc(1)+1)/2 -cm(1);
meshY = meshY - (lenc(2)+1)/2 -cm(2);
meshZ = meshZ - (lenc(3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)

meshX = meshX * ((boxx/h)/lenc(1));
meshY = meshY * ((boxx/h)/lenc(2));
meshZ = meshZ * ((boxx/h)/lenc(3));


rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2); % r^2 cube in Mpc

vrcm=(VcmX.*meshX+VcmY.*meshY+VcmZ.*meshZ)./rcube;

