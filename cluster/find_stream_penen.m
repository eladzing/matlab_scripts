%% plot streams according to conditions
units;

%global SIMTYPE
%global CLUSTER;
global VCM
%global aexpn;
global NCELL;
global zred
global hub;

RVIR=get_rvir();
MVIR=get_mvir();
%VVIR=get_vvir();
TVIR=get_tvir();

vcm = VCM;
cm=[0,0,0];



%% read data

% get velocity 
[hubX,hubY,hubZ] = hubble_flow(boxx,[0,0,0]);
vx = Vx(boxx)+hubX-vcm(1);
vy = Vy(boxx)+hubY-vcm(2);
vz = Vz(boxx)+hubZ-vcm(3);

vc=sqrt(vx.^2+vy.^2+vz.^2);

% get positions
[meshY, meshX, meshZ] = meshgrid(1:size(vy,1), 1:size(vx,2), 1:size(vz,3));
%convert to center origin coordinates
meshX = meshX - (size(vx,1)+1)/2 -cm(1);
meshY = meshY - (size(vy,2)+1)/2 -cm(2);
meshZ = meshZ - (size(vz,3)+1)/2 -cm(3);
% Fix Units (to be in Mpc)
meshX = meshX * ((boxx/hub)/NCELL);
meshY = meshY * ((boxx/hub)/NCELL);
meshZ = meshZ * ((boxx/hub)/NCELL);

rcube=sqrt(meshX.^2+meshY.^2+meshZ.^2) ; % r cube in Mpc

% read density
rog=RHOG(boxx);
[roShell,~]=read_RHO_Profiles(rcube);
ros=RHOTOT(boxx)-RHODM(boxx)-rog;

% read flux
vr=(vx.*meshX+vy.*meshY+vz.*meshZ)./rcube ; %radial velocity

fluxnorm=(0.1117.*(MVIR./1e15).^0.15*(1+zred)^2.25)./(4.*pi); % normalization for flux
[Mg200 , ~, ~]=read_Mass_Profiles(RVIR);
flux=rog.*vr.*rcube.^2./Mg200.*(km/Mpc*Gyr)./fluxnorm ; % normalized



% read temperature
tmp=T(boxx);
tShell=read_T_Profile(rcube);

% mach 
 gm=5/3;
 cso=(gm*kb/mm.*tShell).^0.5;  % ./1e5;  % so
 cso=cso./km;
 
 mach=vc./cso;

 %smooth
smLen=ceil(100./(1e3.*boxx./hub).*NCELL);
smLen=smLen-1+mod(smLen,2);
flux=smooth3(flux,'box',smLen);
mach=smooth3(mach,'box',smLen); %/cso;
tmp=smooth3(tmp,'box',smLen); %/cso;
ross=smooth3(ros,'box',smLen);

% set conditions
fMask=flux<-3;

roMask=rog./roShell>=1;

rosMask=ross<=0;

tmpMask=tmp./tShell<=1;

machMask=mach>1;

%% build radius array

radM=rcube( machMask & fMask & tmpMask);

penen=min(radM)./RVIR