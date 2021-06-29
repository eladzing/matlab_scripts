function [hubX hubY hubZ] = hubble_flow(MPc,cm) 
%%% returns a  cube containing the hubble flow 
%%% velocity in km/sec from the center of the box 
global NCELL;
global hub
global aexp
global zred

if ~exist('cm','var')
    cm = [0,0,0];
end

[meshY, meshX, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);
%convert to center origin coordinates
meshX = meshX - (NCELL+1)/2 -cm(1);
meshY = meshY - (NCELL+1)/2 -cm(2);
meshZ = meshZ - (NCELL+1)/2 -cm(3);
% Fix Units (to be in Mpc)
%h=0.7;
meshX = meshX * ((MPc*aexp/hub)/NCELL);
meshY = meshY * ((MPc*aexp/hub)/NCELL);
meshZ = meshZ * ((MPc*aexp/hub)/NCELL);

Hubble=friedeq(zred);

hubX=Hubble.*meshX;
hubY=Hubble.*meshY;
hubZ=Hubble.*meshZ;
