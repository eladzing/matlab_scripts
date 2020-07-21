% SCALEFRAMES - Scales the frames by their respective weights.
%
% d=scaleFrames(d,W)
% 
% The last dimension of d must be the same as the length of W. W can either
% be an N-by-N matrix or an array of length N.

% By Ran Klein 2009-04-29


% *******************************************************************************************************************
% *                                                                                                                 *
% * Copyright [2014] Ottawa Heart Institute Research Corporation.                                                   *
% * This software is confidential and may not be copied or distributed without the express written consent of the   *
% * Ottawa Heart Institute Research Corporation.                                                                    *
% *                                                                                                                 *
% *******************************************************************************************************************


function d=scaleFrames(d,W)

if ndims(W)==2 && (size(W,1)==1 || size(W,2)==1);
	W = diag(W);
end

s = size(d);
d = reshape(d,[prod(s(1:end-1)) s(end)])*W;
d = reshape(d,s);