function [rvir mvir]=calc_virials(cluster,varargin);


% create density profiles
global NCELL
global hub
global Omm


cellsize=1/hub/Ncell;
ri=1:0.5:(8/(1/NCELL));

rr=ri.*cellsize;


[~ RHOTOT_Profile] = read_RHO_Profiles(rr)

% find critical/mean density of the universe in units of Msun/Mpc^3
%units
rho_critical=2.775e11*(hub)^2;
rho_mean=rho_critical*0mm

rvir1=rvir