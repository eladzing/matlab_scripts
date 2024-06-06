function res = calcCoolingTime(density,energy,coolingRate)
%CALCCoolingTime function to calculate cooling time from internal energy and
% electron abundance in the Illustris simulation. Result is in K.
% Internal Energy is in (km/sec)^2 ( per unit mass)

%units;
global illUnits
global cosmoStruct

eps=(1/cosmoStruct.Hfraction-1)/4;

densCGS=density.*illUnits.densityUnit.*(illUnits.physUnits.Ms/illUnits.physUnits.kpc^3); %density in cgs
energyCGS=energy.*(illUnits.physUnits.km)^2; %internal energy (per mass) in cgs
coolingRateFac=(1+eps)*(1+2*eps).*(cosmoStruct.Hfraction/illUnits.physUnits.mp).^2.*densCGS; 
% simulation cooling rate is lambda/n_H^2. 

res=-1.*energyCGS./(coolingRate.*coolingRateFac)./illUnits.physUnits.Gyr;

res(coolingRate>=0)=0;

end

% 
% 
%   dens_cgs = self.codeDensToPhys(code_dens, cgs=True) # g/cm^3
% 
%         ratefact = self.hydrogen_massfrac**2 / self.mass_proton**2 * dens_cgs # 1/(g*cm^3)
% 
%         coolrate = code_gfmcoolrate * ratefact # erg cm^3/s * (1/g/cm^3) = erg/s/g (i.e. specific rate)
% 
%         u_cgs_spec = code_u * self.UnitVelocity_in_cm_per_s**2 # i.e. (km/s)^2 to (cm/s)^2, so specific erg/g
% 
%         t_cool = u_cgs_spec / (-1.0*coolrate) / self.s_in_Gyr
% 
%  
% 
%         # if lambda_net is positive set t_cool=0 (i.e. actual net heating, perhaps from the background)
% 
%         w = np.where(code_gfmcoolrate >= 0.0)
% 
%         t_cool[w] = 0.0
% 
%  
% 
%         return t_cool
% 
%  


