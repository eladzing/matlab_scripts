
function [VcmX VcmY VcmZ]=Vcm_full(rvcm,cm,hubFactor) 
%returns the components of center of mass velocity calculated within a given radial zone.
%Also takes hubble flow into account

global hub
global aexp
if ~exist('cm')
    cm = [0,0,0];
end
%if ~exist('hubflag')
%    hf=0;
%else
%    hf=(strcmpi(hubflag,'hub'));    
%end

%% calculate Vcm within a given radius

%read velocity and density from appropriate box 
if length(rvcm)<2
    rvcm=[0 rvcm];
end
h=hub;
% find the right box to calculate V_cm
vbox=min(2^max(0,ceil(log2(2*max(rvcm(:))*h/aexp))),8);

[hubX hubY hubZ] = hubble_flow(vbox,cm);

Vxx = Vx(vbox)+hubFactor.*hubX;
Vyy = Vy(vbox)+hubFactor.*hubY;
Vzz = Vz(vbox)+hubFactor.*hubZ;

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



