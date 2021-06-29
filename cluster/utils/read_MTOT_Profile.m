function result = read_MTOT_Profile(R_Profile)
% Utility function to read the cumulative total mass profile
% from Kravtsov's data for the current cluster (See: new_env).
%
% @param R_Profile  The radius profile for which the profile should be
%                   calculated.
% @returns  The total mass profiles
%
error(nargchk(1,1,nargin,'struct'));

global zred
global hub

SKIP_LINES                 = 19;
OUTER_RADIUS_COL           = 3;
TOTAL_CUM_OVER_DENSITY_COL = 10;

% Constants
%Om0 = 0.3;
h = hub;
rhomean=rho_mean(zred)/h^2; % find mean density of the universe
%rhocrit = 2.7755e+11; % [h^2 Msun Mpc^-3].
%rhomean= rhocrit*0.3;

global PROFILE_FILE;
table = dlmread(PROFILE_FILE, '', SKIP_LINES-1, 0);
Outer_Radius_Profile           = table(:, OUTER_RADIUS_COL);
Total_Cum_Over_Density_Profile = table(:, TOTAL_CUM_OVER_DENSITY_COL);

% proper volume enclosed within radius rr [Mpc^3 h^-3]
Volume_Profile = (4.0*pi/3.0)*(Outer_Radius_Profile./(1+zred)).^3;
% total mass within within radius rr [Msun/h]
MTOT_Profile   = (Total_Cum_Over_Density_Profile*rhomean).*Volume_Profile; 


% Interpolate profile to given R_Profile, and convert units to Msun)
result = spline(Outer_Radius_Profile/h/(1+zred), MTOT_Profile/h, R_Profile);