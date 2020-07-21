function result= rho_mean(z,varargin)

%% calculates the mean density of the universe, in units of Msun/Mpc^3
% Calculatopn is done by finding rho_crit*omega_m as a function of redshift

cosmoStructFlag=false;
showFlag=false;

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
        case{'units','show','show units'}
            showFlag=true;
            
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
            error('RHO_MEAN - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


units
if cosmoStructFlag
    rhocrit=3.*(friedeq(z,'cosmoStruct',cosmoStruct).*Units.km./Units.Mpc).^2./(8*pi*Units.G)./Units.Ms*Units.Mpc^3;
    result=rhocrit.*omega_m(z,'cosmoStruct',cosmoStruct);
else
    rhocrit=3.*(friedeq(z,'hub',hub,'omm',omm0,'oml',oml0).*Units.km./Units.Mpc).^2./(8*pi*Units.G)./Units.Ms*Units.Mpc^3;
    result=rhocrit.*omega_m(z,'omm',omm0,'oml',oml0);
end

if showFlag
    fprintf('mean density in units of Msun/Mpc^3 \n');
end

