function res = projectionFac(phi,theta,psi)
%PROJECTIONFAC generate a factor to account for projection effects.
%   Given two roation angles find the projected distance
% 



%phi=pi/4;
%theta=pi/4;

D = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
C = [1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
B = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];

A=B*C*D;

%xVec=[1;0;0];
yVec=[0;1;0];
%zVec=[0;0;1];

%xnew=A*xVec;
ynew=A*yVec;
%znew=A*zVec;

%res(1)=sqrt(sum(xnew(1:2).^2));
res=sqrt(sum(ynew(1:2).^2));
%res(3)=sqrt(sum(znew(1:2).^2));

% Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
% % Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
% Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
% 
% vec=[0;1;0];
% 
% a=Rx*(Rz*vec);
% b=Rz*(Rx*vec);
% res(1)=sqrt(sum(a(1:2).^2));
% res(2)=sqrt(sum(b(1:2).^2));

end

