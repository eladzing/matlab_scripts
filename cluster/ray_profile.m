 function [prof1 prof2 rp]=ray_profile(varargin)

% create profile along a thin line emanating from the center 
% argument pairs:
% boxx - the boxx which to create profile 
% point - vector which defines the direction 
% width - of ray in units of cellsize?
% datatype - which type of data to plot
% data - given cube of data values. 
% weight type - weight of data averaging
% weight cube - given weight datacube


if nargin < 4
    error('ray_profile: not enough arguments')
elseif nargin > 10
    error('ray_profile: too many arguments')
end

width=3;
wtflag=0;
wtypeflag=0;
pointflag=0;
cubeflag=0;
typeflag=0;

RVIR=get_rvir();
MVIR=get_mvir();
%VVIR=get_vvir();
TVIR=get_tvir();


for i=1:2:nargin
    switch varargin{i}
        case {'boxx','box'}
            boxx=varargin{i+1};
        case{'data','datacube'}
            cubeflag=1;
            cube=varargin{i+1};
        case{'point','ray','direction'}
            pointflag=1;
            if(length(point)~=3)
                error('ray_profile: illegal point');
            end
            xp=point(1);yp=point(2);zp=point(3);
        case{'width'}
            width=varargin{i+1};
        case{'weightcube','wtcube','wtdatacube','wt'}
            wtflag=1;
            wt=varargin{i+1};
        case{'datatype','type'}
            typeflag=1;
            type= varargin{i+1};
        case{'weighttype','wttype','wtype'}
            wtypeflag=1;
            wtype= varargin{i+1};
        otherwise
            error('ray_profile: illegal argument %s',varargin{i})
    end
end

if ~any(boxx==[1 2 4 8])
    error('ray_profile: illegal box')
end

if ~(cubeflag || typeflag)
    error('ray_profile: must enter data or datatype')
end
if (cubeflag && typeflag)
    error('ray_profile: too many data  arguments')
end

if typeflag
    switch type 
        case {'T','temp','temperature','Temperature'}
            %cube=T(boxx);
        case {'K','ent','entropy','Entropy'}
            SVIR=TVIR.*RVIR.^2./(3.*MVIR./(4.*pi)).^(2./3); %normalization for entropy
            cube=S(boxx)./SVIR;
        otherwise
                error('ray_profile: unknown data type %s',type);
    end
end

if ~(wtflag || wtypeflag)
    wt=ones(size(cube));
end

if (wtflag && wtypeflag)
    error('ray_profile: too many weight arguments')
end                    

if wtypeflag
    switch wtype 
        case {'rho','mass'}
            wt=RHOG(boxx);
        case {'vol','volume','ones','equal'}
            wt=ones(size(cube));
        otherwise
            error('ray_profile: unknown weight type %s',wtype);
    end
end                
                

global NCELL 
global hub

cellsize=boxx/hub/NCELL;

%% randomly choose a point defining the ray
 
if ~pointflag
    rr=1;
   % select theta & ohi 
   costheta=rand()*2-1;
   phi=rand()*2*pi;
   
   xp=rr*sin(acos(costheta))*cos(phi);
   yp=rr*sin(acos(costheta))*sin(phi);
   zp=rr*costheta;
end

normpoint=sqrt(xp^2+yp^2+zp^2);

%% find all points along the ray

% create coordiante matrix
[meshY, meshX, meshZ] = meshgrid(1:NCELL,1:NCELL,1:NCELL);
        %convert to center origin coordinates
        meshX = meshX - (NCELL+1)/2;
        meshY = meshY - (NCELL+1)/2;
        meshZ = meshZ - (NCELL+1)/2;
        % Fix Units (to be in Mpc)
        meshX = meshX * (boxx/hub/NCELL);
        meshY = meshY * (boxx/hub/NCELL);
        meshZ = meshZ * (boxx/hub/NCELL);

        
%construct t-value matrix        
tval=(meshX.*xp+meshY.*yp+meshZ.*zp)./normpoint^2;

% construct distance matrix 
dist=sqrt((meshY.*zp-meshZ.*yp).^2+(meshX.*zp-meshZ.*xp).^2+(meshX.*yp-meshY.*xp).^2)./normpoint./cellsize;

% construct angle matrix
cos_angl=tval./sqrt(meshX.^2+meshY.^2+meshZ.^2).*normpoint;

clear meshX meshY meshZ

maskp=dist<=width & cos_angl>=0; 
maskm=dist<=width & cos_angl<=0; 

clear dist tval cos_angl

prof1=MAKE_PROFILE_FROM_CUBE_AVG(cube.*maskp,wt.*maskp);
prof2=MAKE_PROFILE_FROM_CUBE_AVG(cube.*maskm,wt.*maskm);
rp=(1:NCELL).*cellsize.*0.5;


