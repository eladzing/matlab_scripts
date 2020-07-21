function  cosmoStruct=set_LCDM_cosmology(varargin)
%SET_LCDM_COSMOLOGY - set generic LCDM model based on Planck2015/2016

verboseFlag=true;

i=1;

while (i<=length(varargin))
    switch lower(varargin{i})
        case{'noshow','quiet'}
            verboseFlag=false;
            
        otherwise
            error('SET_LCDM_COSMOLOGY - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

%global Omm % omega matter
cosmoStruct.Omm = 0.3089;

%global Oml % omega lambda
cosmoStruct.Oml = 0.6911;

%global Omb % omega baryons
cosmoStruct.Omb = 0.0486;

%global hub % hubble parameter in km/sec/MPc / 100
cosmoStruct.hub = 0.6774;

%global sigma8
cosmoStruct.sigma8 = 0.8159;

%global ns % is this the power spectrum slope ?
cosmoStruct.ns = 0.9667;

% hydrogen mass fraction in primordial universe
cosmoStruct.Hfraction = 0.75;


epsilon = (1/cosmoStruct.Hfraction - 1)/4; % helium to hydrogen number density ratio

% mean mass per particle for primordial composition plasma
cosmoStruct.muMass = (1+4* epsilon)/(2+3*epsilon);

% No. of electrons per particle for primordial composition plasma
cosmoStruct.xi = (1+2* epsilon)/(2+3*epsilon);

units;
cosmoStruct.tHubble=1/(100*cosmoStruct.hub*Units.km./Units.Mpc)/Units.Gyr;


if verboseFlag
    fprintf('*** LCDM cosmolgy  has been set: *** \n');
    fprintf('hub= %s \n',num2str(cosmoStruct.hub));
    fprintf('Omm = %s \n',num2str(cosmoStruct.Omm));
    fprintf('Oml= %s \n',num2str(cosmoStruct.Oml));
    fprintf('Omb= %s \n',num2str(cosmoStruct.Omb));
    fprintf('***  *** \n');
end



end

