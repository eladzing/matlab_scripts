function y=rho_nfw(x,ro0,cvir)

% NFW profile. free parameters are ro0 and c_vir
% r must be given in units of R_vir
y=zeros(size(x));
y=ro0./(( x .* cvir) .* (1 + x .* cvir).^2);
end


