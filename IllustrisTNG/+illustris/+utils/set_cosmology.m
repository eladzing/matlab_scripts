function  set_cosmology(mountedFlag )
%SET_COSMOLOGY - set cosmological parameters for the TNG model based on
%Pillepich2017 which in turn is based on Planck2016 results
%

if ~exist('mountedFlag','var')
    mountedFlag=false;
end

global BASEPATH
global cosmoStruct
global simDisplayName

zred0=99;

if contains(simDisplayName,'TNG35')
  zred0=4; 
end

if contains(simDisplayName,'TNG')
    if mountedFlag
        
        
        header=illustris.groupcat.loadHeader(BASEPATH,zred0,0);
        
        %global Omm % omega matter
        cosmoStruct.Omm = header.Omega0; % 0.31;
        
        %global Oml % omega lambda
        cosmoStruct.Oml = header.OmegaLambda;%  0.69;
        
        %global Omb % omega baryons
        cosmoStruct.Omb = 0.0486;
        
        %global hub % hubble parameter in km/sec/MPc / 100
        cosmoStruct.hub = header.HubbleParam;%0.677;
        
        %global sigma8
        cosmoStruct.sigma8 = 0.08159;
        
        %global ns % is this the power spectrum slope ?
        cosmoStruct.ns = 0.97;
        
          % hydrogen mass fraction in primordial universe
        cosmoStruct.Hfraction = 0.76; 
        
    else
        
        %global Omm % omega matter
        cosmoStruct.Omm = 0.31;
        
        %global Oml % omega lambda
        cosmoStruct.Oml = 0.69;
        
        %global Omb % omega baryons
        cosmoStruct.Omb = 0.0486;
        
        %global hub % hubble parameter in km/sec/MPc / 100
        cosmoStruct.hub = 0.677;
        
        %global sigma8
        cosmoStruct.sigma8 = 0.08159;
        
        %global ns % is this the power spectrum slope ?
        cosmoStruct.ns = 0.97;
        
        % hydrogen mass fraction in primordial universe
        cosmoStruct.Hfraction = 0.76; 
    end
    
elseif contains(simDisplayName,'Illustris')
    %global Omm % omega matter
    cosmoStruct.Omm = 0.2726;
    
    %global Oml % omega lambda
    cosmoStruct.Oml = 0.7274;
    
    %global Omb % omega baryons
    cosmoStruct.Omb = 0.0456;
    
    %global hub % hubble parameter in km/sec/MPc / 100
    cosmoStruct.hub = 0.704;
    
      
    % hydrogen mass fraction in primordial universe
        cosmoStruct.Hfraction = 0.76; 
        
    %     %global sigma8
    %     cosmoStruct.sigma8 = 0.08159;
    %
    %     %global ns % is this the power spectrum slope ?
    %     cosmoStruct.ns = 0.97;
    
else
    error('SET_COSMOLOGY - Unknown Simulation: %s',simDisplayName)
end



epsilon = (1/cosmoStruct.Hfraction - 1)/4; % helium to hydrogen number density ratio

% mean mass per particle for primordial composition plasma
cosmoStruct.muMass = (1+4* epsilon)/(2+3*epsilon);  

% No. of electrons per particle for primordial composition plasma
cosmoStruct.xi = (1+2* epsilon)/(2+3*epsilon);  

units;
cosmoStruct.tHubble=1/(100*cosmoStruct.hub*Units.km./Units.Mpc)/Units.Gyr;


fprintf('*** cosmolgy for %s has been set: *** \n',simDisplayName);
fprintf('hub= %s \n',num2str(cosmoStruct.hub));
fprintf('Omm = %s \n',num2str(cosmoStruct.Omm));
fprintf('Oml= %s \n',num2str(cosmoStruct.Oml));
fprintf('Omb= %s \n',num2str(cosmoStruct.Omb));
fprintf('***  *** \n');




end

