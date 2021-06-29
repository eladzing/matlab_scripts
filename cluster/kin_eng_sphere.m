function [Ekin_sph Ekin_r_sph]= kin_eng_sphere(boxx,vcmflag,rvcm,hubflag) 

%Calculates the kinetic energy of the gas and maps it on a spherical cube
%Does the same for the kinetuc energy associated witht the radial velocity

[Ek Ekr]=kin_eng(boxx,vcmflag,rvcm,hubflag);

Ekin_sph = cart2sphere(Ek);
Ekin_r_sph=cart2sphere(Ekr);
