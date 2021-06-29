function vrcm_sph=Vrcm_SPHERE(boxx,rvcm,cm,hubflag) 
%returns a spherical cube of the radial component of
%Vcm calculated within the mask 'ind'
if ~exist('cm')
    cm = [0,0,0];
end
vrcm_sph=cart2sphere(Vrcm(boxx,rvcm,cm,hubflag));



