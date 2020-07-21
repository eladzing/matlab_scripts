function div=calc_divV_single(Vxx,Vyy,Vzz,boxx)

% calculate the divergence of the cartesian velocity field
% assumes velocities are already 'fixed', i.e. with hubble flow
% and vcm corrected.

global NCELL
global hub

[meshX, meshY, meshZ] = meshgrid(1:NCELL, 1:NCELL, 1:NCELL);

%convert to center origin coordinates
        meshX = meshX - (NCELL+1)/2;% -cm(1);
        meshY = meshY - (NCELL+1)/2;% -cm(2);
        meshZ = meshZ - (NCELL+1)/2;% -cm(3);
        % Fix Units (to be in Mpc)
        meshX = meshX * ((boxx/hub)/NCELL);
        meshY = meshY * ((boxx/hub)/NCELL);
        meshZ = meshZ * ((boxx/hub)/NCELL);
    
        div=divergence_single(meshX,meshY,meshZ,Vxx,Vyy,Vzz);
end

function div=divergence_single(varargin)
%DIVERGENCE  Divergence of a vector field.
%   DIV = DIVERGENCE(X,Y,Z,U,V,W) computes the divergence of a 3-D
%   vector field U,V,W. The arrays X,Y,Z define the coordinates for
%   U,V,W and must be monotonic and 3-D plaid (as if produced by
%   MESHGRID).
%   
%   DIV = DIVERGENCE(U,V,W) assumes 
%         [X Y Z] = meshgrid(1:N, 1:M, 1:P) where [M,N,P]=SIZE(U). 
%
%   DIV = DIVERGENCE(X,Y,U,V) computes the divergence of a 2-D
%   vector field U,V. The arrays X,Y define the coordinates for U,V
%   and must be monotonic and 2-D plaid (as if produced by
%   MESHGRID). 
%   
%   DIV = DIVERGENCE(U,V) assumes 
%         [X Y] = meshgrid(1:N, 1:M) where [M,N]=SIZE(U). 
%   
%   Example:
%      load wind
%      div = divergence(x,y,z,u,v,w);
%      slice(x,y,z,div,[90 134],[59],[0]); shading interp
%      daspect([1 1 1])
%      camlight
%
%   See also STREAMTUBE, CURL, ISOSURFACE.

%   Copyright 1984-2009 The MathWorks, Inc. 
%   $Revision: 1.4.4.3 $  $Date: 2011/03/09 07:03:59 $

error(nargchk(2,6,nargin,'struct'));
[x y z u v w] = parseargs(nargin,varargin);

% Take this out when other data types are handled
%u = double(u);
%v = double(v);
%w = double(w);

if isempty(w)  % 2D

  [msg x y] = xyuvcheck_local(x,y,u,v);  error(msg) 
  if isempty(x)
    [px junk] = gradient(u); %#ok
    [junk qy] = gradient(v); %#ok
  else
    hx = x(1,:); 
    hy = y(:,1); 
    [px junk] = gradient(u, hx, hy); %#ok
    [junk qy] = gradient(v, hx, hy); %#ok
  end
  div = px+qy;
  
else   %3D
  
  [msg x y z] = xyzuvwcheck_local(x,y,z,u,v,w);  error(msg) 
  if isempty(x)
    [px junk junk] = gradient(u); %#ok
    [junk qy junk] = gradient(v); %#ok
    [junk junk rz] = gradient(w); %#ok
  else
    hx = x(1,:,1); 
    hy = y(:,1,1); 
    hz = z(1,1,:); 
    [px junk junk] = gradient(u, hx, hy, hz); %#ok
    [junk qy junk] = gradient(v, hx, hy, hz); %#ok
    [junk junk rz] = gradient(w, hx, hy, hz); %#ok
  end
  
  div = px+qy+rz;
  
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, z, u, v, w] = parseargs(nin, vargin)

x = [];
y = [];
z = [];
w = [];

if nin==2         % divergence(u,v)
  u = vargin{1};
  v = vargin{2};
elseif nin==3     % divergence(u,v,w)
  u = vargin{1};
  v = vargin{2};
  w = vargin{3};
elseif nin==4     % divergence(x,y,u,v)
  x = vargin{1};
  y = vargin{2};
  u = vargin{3};
  v = vargin{4};
elseif nin==6     % divergence(x,y,z,u,v,w)
  x = vargin{1};
  y = vargin{2};
  z = vargin{3};
  u = vargin{4};
  v = vargin{5};
  w = vargin{6};
else
  error(message('MATLAB:divergence:WrongNumberOfInputs')); 
end
end


function [msg,nx,ny,nz] = xyzuvwcheck_local(x,y,z,u,v,w)
%XYZUVWcheck_local  Check arguments to 3D vector data routines.
%   [MSG,X,Y,Z] = XYZUVWCHECK(X,Y,Z,U,V,W) checks the input arguments
%   and returns either an error message structure in MSG or valid 
%   X,Y,Z. The ERROR function describes the format and use of the
%   error structure.
%
%   See also ERROR

%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 1.6.4.2 $  $Date: 2009/11/13 04:37:57 $

msg = struct([]);
nx = x;
ny = y;
nz = z;

sz = size(u);
if ~isequal(size(v), sz) || ~isequal(size(w), sz)
  msg = makemsg('UVWSizeMismatch','U,V,W must all be the same size.');
  return
end

if ndims(u)~=3
  msg = makemsg('UVWNot3D','U,V,W must all be a 3D arrays.');
  return
end
if min(sz)<2
  msg = makemsg('UVWPlanar','U,V,W must all be size 2x2x2 or greater.'); 
  return
end

nonempty = ~[isempty(x) isempty(y) isempty(z)];
if any(nonempty) && ~all(nonempty)
  msg = makemsg('XYZMixedEmpty','X,Y,Z must all be empty or all non-empty.');
  return;
end

if ~isempty(nx) && ~isequal(size(nx), sz)
  nx = nx(:);
  if length(nx)~=sz(2)
    msg = makemsg('XUSizeMismatch','The size of X must match the size of U or the number of columns of U.');
    return
  else
    nx = repmat(nx',[sz(1) 1 sz(3)]);
  end
end

if ~isempty(ny) && ~isequal(size(ny), sz)
  ny = ny(:);
  if length(ny)~=sz(1)
    msg = makemsg('YUSizeMismatch','The size of Y must match the size of U or the number of rows of U.');
    return
  else
    ny = repmat(ny,[1 sz(2) sz(3)]);
  end
end

if ~isempty(nz) && ~isequal(size(nz), sz)
  nz = nz(:);
  if length(nz)~=sz(3)
    msg = makemsg('ZUSizeMismatch','The size of Z must match the size of U or the number of pages of U.');
    return
  else
    nz = repmat(reshape(nz,[1 1 length(nz)]),[sz(1) sz(2) 1]);
  end
end
end

function msg = makemsg(id,str)
msg.identifier = ['MATLAB:xyzuvwcheck:' id];
msg.message = str;
end