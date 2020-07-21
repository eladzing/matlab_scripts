function full_flux=new_flux(boxx,varargin)
 
%% despite a bunch of other files dealing with flux -I decided to write 
% this one rather than try and remember who does what 
%
%This file calculatees the full mass flux for a given box, 
%
%in units of  Msun/yr

%%calculate radial velocity and transform into sphere
units 
vrflag=false;
roflag=false;
i=1;
while i<=nargin-1
    
    switch varargin{i}
        case 'vr'
            i=i+1;
            vr=varargin{i};
            vr_sph=cart2sphere(vr);
            vrflag=true;
        case 'vr_sph'
             i=i+1;
             vr_sph=varargin{i};
             vrflag=true;
        case 'ro_sph'
             i=i+1;
             ro_sph=varargin{i};
             roflag=true;
        otherwise
            error('new_flux: Illegal  option')
    end
    i=i+1;
end

if ~vrflag 
    vr_sph=cart2sphere(Vr_full(boxx));
end
if ~roflag
    ro_sph=RHOG_sphere(boxx);
end


ds=ds_sphere(boxx);

%% factor calculation: 
% [ro]=Msun/Mpc^3 
%  [vr]=km/sec
%  [ds]=Mpc^2

fac=km/Mpc*yr; 

full_flux=ro_sph.*vr_sph.*ds.*fac;

