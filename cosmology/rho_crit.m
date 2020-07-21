function result= rho_crit(z,varargin)

%% calculates the mean density of the universe, in units of Msun/Mpc^3 
% Calculatopn is done by finding rho_crit as a function of redshift


cosmoStructFlag=false;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'hub','hubble','h'}
            i=i+1;
            hub=varargin{i};
        case{'omm','o_m','omegam','omega_m'}
            i=i+1;
            omm0=varargin{i};
        case{'oml','o_l','omegal','omega_l'}
            i=i+1;
            oml0=varargin{i};
        case{'cosmostruct','cosmo','cosmology'}
            i=i+1;
            cosmoStruct=varargin{i};
            if ~isstruct(cosmoStruct)
                error('OMEGA_M: cosmology structure must be a structure')
            end
            cosmoStructFlag=true;
            %                omm0=cosmoStruct.omm;
            %                oml0=cosmoStruct.oml;
            %                hub=cosmoStruct.hub;
        otherwise
            error('OMEGA_M - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


units 

if cosmoStructFlag
    result=3.*(friedeq(z,'cosmoStruct',cosmoStruct).*Units.km./Units.Mpc).^2./(8*pi*Units.G)./Units.Ms*Units.Mpc^3;
   
else
    result=3.*(friedeq(z,'hub',hub,'omm',omm0,'oml',omm0).*Units.km./Units.Mpc).^2./(8*pi*Units.G)./Units.Ms*Units.Mpc^3;
   
end

