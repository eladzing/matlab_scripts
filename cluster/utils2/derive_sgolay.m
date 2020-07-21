function [derivative n_range] = derive_sgolay(data, dx, ORDER, WINDOW_LENGTH)
% Savitzky-Golay first derivative, works on spherical 3D data (derives along R), and vectors.
% See: derive_simple().

[b,g] = sgolay(ORDER,WINDOW_LENGTH);   % Calculate S-G coefficients

HalfWin  = ((WINDOW_LENGTH+1)/2) -1;
n_range = (WINDOW_LENGTH+1)/2:length(data)-(WINDOW_LENGTH+1)/2;

SG2 = zeros([max(n_range), size(data, 2), size(data, 3)]);
for n = n_range
%     for m = 1:size(data,2)
%         for p = 1:size(data,3)
%             % 1st differential
%             SG1(n,m,p) = dot(g(:,2), data(n - HalfWin: n + HalfWin,m,p));
%         end
%     end
    if (ndims(data) == 3)
        SG2(n,:,:) = dot(repmat(g(:,2), [1 size(data,2) size(data,3)]), data(n - HalfWin: n + HalfWin, :, :), 1);
    elseif (isvector(data))
        SG2(n) = dot(g(:,2), data(n - HalfWin: n + HalfWin));
    else
        error('Wrong number of dimentios!');
    end
end


% Turn differential into derivative
if (ndims(data) == 3)
    derivative = SG2(n_range,:,:)/dx;
elseif (isvector(data))
    derivative = SG2(n_range)/dx;
else
    error('Wrong number of dimentios!');
end

