function [RHOG_Profile,RHOTOT_Profile] = read_RHO_Profiles(R_Profile)
% Utility function to read the gas and total density profiles
% from Kravtsov's data for the current cluster (See: new_env).
%
% @param R_Profile  The radius profile for which the profiles should be
%                   calculated, in proper units 
% @returns  The gas, and total density profiles
%
narginchk(1,1);

global zred
global hub

SKIP_LINES            = 19;
INNER_RADIUS_COL      = 2;
OUTER_RADIUS_COL      = 3;
DIFF_GAS_DENSITY_COL  = 4;
DIFF_STAR_DENSITY_COL = 5;
DIFF_DM_DENSITY_COL   = 6;

% Constants
%Om0 = 0.3;
h = hub;
rhomean=rho_mean(zred)/h^2; % find mean density of the universe
%rhocrit = 2.7755e+11; % [h^2 Msun Mpc^-3].

global PROFILE_FILE;
table = dlmread(PROFILE_FILE, '', SKIP_LINES-1, 0);
Inner_Radius_Profile       = table(:, INNER_RADIUS_COL);
Outer_Radius_Profile       = table(:, OUTER_RADIUS_COL);
Diff_Gas_Density_Profile   = table(:, DIFF_GAS_DENSITY_COL);
Diff_Total_Density_Profile = sum(table(:, [DIFF_GAS_DENSITY_COL, DIFF_STAR_DENSITY_COL, DIFF_DM_DENSITY_COL]), 2);

% Interpolate profile to given R_Profile, and convert units to Msun.
RHOG_Profile = spline((Inner_Radius_Profile + Outer_Radius_Profile)/h/(1+zred)/2, (Diff_Gas_Density_Profile*rhomean)*h^2, R_Profile);
RHOTOT_Profile = spline((Inner_Radius_Profile + Outer_Radius_Profile)/h/(1+zred)/2, (Diff_Total_Density_Profile*rhomean)*h^2, R_Profile);