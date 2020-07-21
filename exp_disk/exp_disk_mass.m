function res=exp_disk_mass(r,beta)
%% unitless part of exponential disk mass profile 
% to get the true profile one needs to multiply by the total disk mass
% r is in units of rd
if ~exist('beta','var')
    beta=1;
end

res=1-exp(-1.*r.*beta).*(1+r.*beta);
end