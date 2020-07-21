function result = lambda_Profile(R_Profile, J_Profile)
% Calculates the scalar spin parameter 
% \lambda = (J/M_g)/(rV_c) where V_c^2=GM_tot(r)/r
% for the current cluster (See: new_env).
%
% @param R_Profile  The radius profile matching the given J_Profile.
% @param J_Profile  The Angular Momentum profile (See:
%                   AngularMomentum_Profiles).
error(nargchk(2,2,nargin,'struct'));

units;

global HALO_PATH
MTOT_Profile = read_MTOT_Profile(HALO_PATH, R_Profile);
[MG_Profile MSTAR_Profile MDM_Profile] = read_Mass_Profiles(HALO_PATH, R_Profile);
Vc = sqrt(G * (MTOT_Profile * Msun) ./ (R_Profile * MPc));
result = ((J_Profile *(1e3*MPc*Msun)) ./ (MG_Profile * Msun)) ./ ((R_Profile * MPc) .* Vc);
