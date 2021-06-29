function [MG_Profile, MSTAR_Profile, MDM_Profile] = read_Mass_Profiles(R_Profile)
% Utility function to read the cumulative gas, star, and DM mass profiles
% from Kravtsov's data for the current cluster (See: new_env).
%
% @param R_Profile  The radius profile for which the profiles should be
%                   calculated.
% @returns  The gas, star, and DM mass profiles
%
narginchk(1,1);

global zred
global hub
global Oml
global Omm

SKIP_LINES                 = 19;
OUTER_RADIUS_COL           = 3;
GAS_CUM_OVER_DENSITY_COL   = 7;
STAR_CUM_OVER_DENSITY_COL  = 8;
DM_CUM_OVER_DENSITY_COL    = 9;

% Constants
%Om0 = 0.3;
h = hub;
rhomean=rho_mean(zred,'hub',hub,'Omm',Omm,'Oml',Oml)/h^2; % find mean density of the universe
%rhocrit = 2.7755e+11; % [h^2 Msun Mpc^-3].


global PROFILE_FILE;
table = dlmread(PROFILE_FILE, '', SKIP_LINES-1, 0);
Outer_Radius_Profile           = table(:, OUTER_RADIUS_COL);
Gas_Cum_Over_Density_Profile   = table(:, GAS_CUM_OVER_DENSITY_COL);
Star_Cum_Over_Density_Profile  = table(:, STAR_CUM_OVER_DENSITY_COL);
DM_Cum_Over_Density_Profile    = table(:, DM_CUM_OVER_DENSITY_COL);

% proper volume enclosed within radius rr [Mpc^3 h^-3]
Volume_Profile = (4.0*pi/3.0)*(Outer_Radius_Profile./(1+zred)).^3;
% total mass within within radius rr [Msun/h]
MG_Profile    = (Gas_Cum_Over_Density_Profile*rhomean).*Volume_Profile; 
MSTAR_Profile = (Star_Cum_Over_Density_Profile*rhomean).*Volume_Profile; 
MDM_Profile   = (DM_Cum_Over_Density_Profile*rhomean).*Volume_Profile; 

% Interpolate profile to given R_Profile, and convert units to Msun)
MG_Profile    = spline(Outer_Radius_Profile/h/(1+zred), MG_Profile/h, R_Profile);
MSTAR_Profile = spline(Outer_Radius_Profile/h/(1+zred), MSTAR_Profile/h, R_Profile);
MDM_Profile   = spline(Outer_Radius_Profile/h/(1+zred), MDM_Profile/h, R_Profile);
