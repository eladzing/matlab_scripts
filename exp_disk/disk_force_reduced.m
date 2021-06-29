function res=disk_force_reduced(r,varargin)
%%function for calculating the unitless part of the gravtiational force
% per unit area acting on the gas in an exponential disk
% r    is in units of the scale radius rd unless stated otherwisw
% beta is the ratio between rd and rg (star & gas scale radii), for a case
%      of different exponential distributions of the gas and stars
% fg   is the ratio of gas to stellar mass.
% fb   is the ratio of bulge to stellar mass


rd=1;
beta=1;
fg=0.1;
fb=0;
BT=0;
md=Inf; %stellar to halo mass ratio (includes bulge) defualt - no halo 
cv=10;
rs=0.01;

starflag=false;% neglects the effect of the gas in determining the potential
bulge='hern'; % bulge model
halo='nfw';
xi=100;  % bulge scale radius

i=1;
while i<=length(varargin)
    switch(varargin{i})
        case 'beta'
            i=i+1;
            beta=varargin{i};
        case {'fg'}
            i=i+1;
            fg=varargin{i};
        case {'fb'}
            i=i+1;
            fb=varargin{i};
        case {'BT' 'B/T'}
            i=i+1;
            BT=varargin{i};
            fb=BT/(1-BT);
        case {'rd'}
            i=i+1;
            rd=varargin{i};
        case {'xi'}
            i=i+1;
            xi=varargin{i};
        case {'bulge'}
            i=i+1;
            bulge=varargin{i};
        case {'star','stellar'}
            starflag=true;
          case {'halo'}
            i=i+1;
            bulge=varargin{i};  
            case {'md'}
            i=i+1;
            bulge=varargin{i};
            case {'cv'}
            i=i+1;
            bulge=varargin{i};
        otherwise
            error('disk_force_reduced: Illegal option')
    end
    i=i+1;
end

if starflag
    fg=0;
end
r=r./rd;

mvms=(1+fb)./md;

%res=beta.^2.* exp(-1.*r*beta).*(expdisk_accel(r,'fg',fg,'beta',beta)+ 2*fb.*bulge_accel(r,'point'));

gDisk=expdisk_accel(r,'fg',fg,'beta',beta);
gBulge=bulge_accel(r,bulge,'fb',fb,'xi',xi);
gHalo=halo_accel(r,mvms,rs,halo,'cv',cv);

res=expdisk_density(r,'gas','beta',beta,'fg',1.0).*...
    (gDisk+gBulge+gHalo);

