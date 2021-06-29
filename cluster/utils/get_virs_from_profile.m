function [rvir, mvir] = get_virs_from_profile(delta_vir)
% Utility function to find the virial radius given a required overdensity
% from Kravtsov's data for the current cluster (See: new_env).
%
% @param delta_vir  The overdensity for which we wih to find Rvir
%                 
% @returns  The radius in Mpc proper. 
%
narginchk(1,1);

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


global PROFILE_FILE;
table = dlmread(PROFILE_FILE, '', SKIP_LINES-1, 0);
Outer_Radius_Profile           = table(:, OUTER_RADIUS_COL);
Total_Cum_Over_Density_Profile = table(:, TOTAL_CUM_OVER_DENSITY_COL);

%interpolate the radius
ind=find(Total_Cum_Over_Density_Profile<1e5);
rv=spline(Total_Cum_Over_Density_Profile(ind),Outer_Radius_Profile(ind),delta_vir);

% proper volume enclosed within radius rr [Mpc^3 h^-3]
Volume_Profile = (4.0*pi/3.0)*(Outer_Radius_Profile./(1+zred)).^3;
% % comoving volume enclosed within radius rr [Mpc^3 h^-3]
% Volume_Profile = (4.0*pi/3.0)*Outer_Radius_Profile.^3;

% total mass within within radius rr [Msun/h]
MTOT_Profile   = (Total_Cum_Over_Density_Profile*rhomean).*Volume_Profile; 

% Interpolate profile to given R_Profile, and convert units to Msun)
mvir=spline(Outer_Radius_Profile, MTOT_Profile/h,rv);

rvir=rv./h/(1+zred); 
