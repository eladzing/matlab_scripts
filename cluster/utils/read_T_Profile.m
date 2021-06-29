function result = read_T_Profile(R_Profile)
% Utility function to read the temperature profile
% from Kravtsov's data for the current cluster (See: new_env).
%
% @param R_Profile  The radius profile for which the profile should be
%                   calculated.
% @returns  The temperature profiles
%
narginchk(1,1);

SKIP_LINES                 = 19;
MASS_WEIGHTED_RADIUS_COL   = 1;
MASS_WEIGHTED_TEMP_COL     = 11;
h = 0.7;

global PROFILE_FILE;
table = dlmread(PROFILE_FILE, '', SKIP_LINES-1, 0);
Mass_Weighted_Radius_Profile = table(:, MASS_WEIGHTED_RADIUS_COL);
Mass_Weighted_Temp_Profile   = table(:, MASS_WEIGHTED_TEMP_COL);

global zred

% Interpolate profile to given R_Profile, and convert units to Msun)
result = spline(Mass_Weighted_Radius_Profile/h/(1+zred), Mass_Weighted_Temp_Profile, R_Profile);