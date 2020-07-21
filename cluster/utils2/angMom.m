function [Jx, Jy, Jz]=angMom(boxx,varargin)
 
%% _function to find the specific angular momentum of a cell
% units are in Mpc * km /sec

global NCELL 
global hub
center = [0,0,0];

i=1;
while i<=length(varargin)
    switch varargin{i}
        case 'center'
            i=i+1;
            center=varargin{i};
        otherwise
            error('angMom: Illegal argument %s',varargin{i});
    end
    i=i+1;
end



vv=zeros(NCELL,NCELL,NCELL,3);
rr=vv;


%% get velocities and create 4D array 
[vx, vy, vz] = get_velocities(boxx);
vv(:,:,:,1)=vx; 
vv(:,:,:,2)=vy;
vv(:,:,:,3)=vz;
%clear vx vy vz;
%% get positions and create 4D array 

[meshY, meshX, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);

%convert to center origin coordinates
meshX = meshX - (NCELL+1)/2 -center(1);
meshY = meshY - (NCELL+1)/2 -center(2);
meshZ = meshZ - (NCELL+1)/2 -center(3);
% Fix Units (to be in Mpc)
meshX = meshX * ((boxx/hub)/NCELL);
meshY = meshY * ((boxx/hub)/NCELL);
meshZ = meshZ * ((boxx/hub)/NCELL);
 
rr(:,:,:,1)=meshX;
rr(:,:,:,2)=meshY;
rr(:,:,:,3)=meshZ;
%clear meshX meshY meshZ;

jj=cross(vv,rr);

Jx=jj(:,:,:,1);
Jy=jj(:,:,:,2);
Jz=jj(:,:,:,3);

%% test
jx2=(vy.*meshZ-vz.*meshY);
jy2=(vz.*meshX-vx.*meshZ);
jz2=(vx.*meshY-vy.*meshX);

i=1

% 
% %%calculate radial velocity and transform into sphere
% %units 
% %vrflag=false;
% %roflag=false;
% i=1;
% while i<=nargin-1
%     
%     switch varargin{i}
%         case 'vr'
%             i=i+1;
%             vr=varargin{i};
%             vr_sph=cart2sphere(vr);
%             vrflag=true;
%         case 'vr_sph'
%              i=i+1;
%              vr_sph=varargin{i};
%              vrflag=true;
%         case 'ro_sph'
%              i=i+1;
%              ro_sph=varargin{i};
%              roflag=true;
%         otherwise
%             error('new_flux: Illegal  option')
%     end
%     i=i+1;
% end
% 
% %if ~vrflag 
% %    vr_sph=cart2sphere(Vr_full(boxx));
% %end
% %if ~roflag
% %    ro_sph=RHOG_sphere(boxx);
% %end
% 
% 
% %ds=ds_sphere(boxx);
% 
% %% factor calculation: 
% % [ro]=Msun/Mpc^3 
% %  [vr]=km/sec
% %  [ds]=Mpc^2
% 
% fac=km/Mpc*yr; 
% 
% full_flux=ro_sph.*vr_sph.*ds.*fac;
% 
