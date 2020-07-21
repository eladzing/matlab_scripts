function [Vxx, Vyy, Vzz] = get_velocities(boxx,varargin) %,vcmflag,rvcm,hubflag,center)
% Returns the velocity field after correcting for center of mass velocity
% and hubble flow.
% Vcm is either calculater

global VCM

%defaults
center = [0,0,0];
hubFactor=true;  %hubble flag
rvcm=[0 1];
findVcmFlag=false;


i=1;
while i<=length(varargin)
    switch varargin{i}
        case 'center'
            i=i+1;
            center=varargin{i};
        case 'nohub'
            hubFactor=0;
        case 'rvcm' 
            i=i+1;
            rvcm=varargin{i};
        case 'findVCM'
            findVcmFlag=true;
        otherwise
            error('get_velocities: Illegal argument %s',varargin{i});
    end
    i=i+1;
end
   
if findVcmFlag
    [VcmX, VcmY, VcmZ]=Vcm_full(rvcm,center,hubFactor);
else
    VcmX=VCM(1); VcmY=VCM(2); VcmZ=VCM(3); % zeros(size(rog)); VcmY=zeros(size(rog)); VcmZ=zeros(size(rog));
end
clear rvc;

[hubX, hubY, hubZ] = hubble_flow(boxx,center);
Vxx = Vx(boxx)+hubFactor.*hubX-VcmX;
Vyy = Vy(boxx)+hubFactor.*hubY-VcmY;
Vzz = Vz(boxx)+hubFactor.*hubZ-VcmZ;


