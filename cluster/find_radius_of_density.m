function res=find_radius_of_density(rhoValue,varargin)
%% find the radius which demarks a given density value (in regular or number denstiy)
% result is in Mpc 

units;


fac=1; 
i=1;
while i<=length(varargin)
    switch lower(varargin{i});
        case {'number','n'}
            fac=Ms/mm/Mpc^3;
        otherwise
            error('find_radius_of_density: Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

r=10.^(-2:0.01:0.5).*get_rvir;
[rop,~] = read_RHO_Profiles(r);

rop=rop.*fac;

res=interp1(rop,r,rhoValue);

end





