function result= omega_m(z,varargin)

%% calculates Omega matter as a function of redshift.
%optional arguments for designating omega_m and omega_l at z=0

%% default
omm0=0.3;
oml0=0.7;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
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
            omm0=cosmoStruct.Omm;
            oml0=cosmoStruct.Oml;
            
        otherwise
            error('OMEGA_M - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end

result= omm0.*(1+z).^3./(oml0+omm0.*(1+z).^3);


