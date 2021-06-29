function result = xrayProj(proj,boxx,zemit,varargin)
% Utility function to load the temperature cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256 x 256 x ne matrix of singles.
% ne is the number of energy bins in the spectrum .
% result is in units of counts / sec / keV
%
%narginchk(3,4);
baseName='xrayFlux_%s_z%s';


i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case {'co','cutout'}
            baseName='xrayFlux_%s_z%s_CO';
        case {'reg','regular','full'}
            baseName='xrayFlux_%s_z%s';
        otherwise
            error('xrayProj - Illegal argument: %s',varargin{i})
    end
    i=i+1;
end


flipFlag=true;
proj0=proj;

%% deal with projection mismatch
switch(lower(proj))
    case{'xy','yx'}
        proj='yz';
        
    case{'yz','zy'}
        proj='xy';
        flipFlag=false;
        
    case{'xz','zx'}
        proj='xz';
        
    otherwise
        error('xrayProj - Illegal projection: %s',proj)
end
global FILE_FORMAT_XRAYPROJ;



result = read_xray_proj(sprintf(FILE_FORMAT_XRAYPROJ,sprintf(baseName,lower(proj),num2str(zemit)), boxx));
%%result = read_xray_proj(sprintf(FILE_FORMAT_XRAYPROJ,sprintf('xrayFlux_%s_z%s_CO',lower(proj),num2str(zemit)), boxx));

if flipFlag
    result.data=permute(result.data,[3 2 1]);
else
    result.data=permute(result.data,[2 3 1]);
end
result.projPlane=proj0;
result.box=boxx;
result.zEmit=zemit; %redshift of emmiting source.
%% correct for volume
% global NCELL
% global hub
% units;
% cellVol=(boxx/hub/NCELL.*Mpc).^3;
% result.data2=result.data.*cellVol;

end
%result = result.data;
