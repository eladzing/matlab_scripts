function res=disk_force_reducede(r,beta)
%% auxilary function to calculate the 'bessel' part of the  
% gravitational acceleration of an exponential disk.
% r is in units of the scale radius rd 

if ~exist('beta','var')
    beta=1;
end

res= exp(-1.*r/beta).*r.*(besseli(0,0.5.*r).*besselk(0,0.5.*r)-besseli(1,0.5.*r).*besselk(1,0.5.*r));

end 