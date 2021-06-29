function result = xrayTempProj(proj,boxx)
% Utility function to load the emmision weighted temperature cube
%
% @param MPSec The required Mpc resolution (1, 2, 4, or 8)
% @returns  A 256 x 256 x ne matrix of singles. 
% ne is the number of energy bins in the spectrum . 
% result is in units of counts / sec / keV 
%
flipFlag=true;
proj0=proj;
narginchk(2,2);
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

result = read_xray_proj_temp(sprintf(FILE_FORMAT_XRAYPROJ,sprintf('xrayTemp_%s',lower(proj)), boxx));


if flipFlag
    result.proj=permute(result.proj,[2 1]);
   % dataW=permute(resultW.data,[3 2 1]);
end
%    dataT=permute(resultT.data,[2 3]);
    %dataW=permute(resultW.data,[2 3 1]);
%end

%resultT.proj=dataT
result.projPlane=proj0;





end
%result = result.data;
