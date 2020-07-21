function res = partialSphereVolume(r1,r2,Radius,sphereCenter)
%PARTIALSPHEREVOLUME - calculates the fractional volume of a sphere
%enclosed between the two radii r1 and r2
%   given a sphere ith center sphereCenter and radius Radius, find the
%   fraction of the sphere volume which is enclosed between two  parallel
%   planes found at r1 and r2


% if size(r1)~=size(r2)
%     error('PARTIALSPHEREVOLUME - integration limits r1 r2 must be of equal length')
% end
% 
% 
% if ~exist('sphereCenter','var')
%     sphereCenter=0;
% end

% % arrange range so r1 is less than r2
% if r1>r2
%     rt=r1;
%     r1=r2;
%     r2=rt;
% end
% 

r1=double(r1./Radius);
r2=double(r2./Radius);
rc=double(sphereCenter./Radius);

r1=max(r1,rc-1);
r2=min(r2,rc+1);

res=0.75.*( (1-rc.^2).*(r2-r1)-1/3.*(r2.^3-r1.^3)+rc.*(r2.^2-r1.^2));

end

