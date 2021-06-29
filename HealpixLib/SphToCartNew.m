% Spherical-Cartesian Coordinate Conversion

function OUT = SphToCartNew(IN,rho)

if ~exist('rho','var')
    rho=1;
end

n = size(IN, 1);

OUT = zeros(n, 3);

phi = IN(:, 1);
theta = IN(:, 2);
r=rho.*sin(phi);

OUT(:,1)=r.*cos(theta);
OUT(:,2)=r.*sin(theta);
OUT(:,3)=rho.*cos(phi);

%for i = 1:n
%    phi = IN(i, 1);
%    theta = IN(i, 2);
    %rho = 1;
    
    %r = rho * sin(phi);
    %OUT(i, :) = [r * cos(theta), r * sin(theta), rho * cos(phi)];
%end

end