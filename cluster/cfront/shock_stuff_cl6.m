%% play around wih shocks 

units;

%for cl6 
% shocks:left to right 
prat=[5.275/2.942 4.366/1.708 1.54/0.5624  ];
u=[1.3 0.1041 -0.4338]*get_vvir;
shkSign=[-1 1 1];
tm=[0.7136 1.495 0.9468].*get_tvir; 
cfT=[0.9112./2.153] 
cfRo=[5.401/2.268]

aa=cfT.*cfRo
%% 1st shock 
%find sonic velocity 
gamma=5/3;
cs=sqrt(gamma*kb/mu/mp.*tm)./1e5;  % sonic velocity in km/sec 

% find shock velocity 


%ushock1 = u1+cs*sqrt(1+(gamma+1)/(2*gamma)*(prat-1))
mach1= sqrt((prat*(gamma+1)+(gamma-1))/(2*gamma))
ushock=u+shkSign.*cs.*mach1

% 
% %% 2nd shock 
% % find sonic velocity 
% gamma=5/3;
% tm=0.985;
% cs=sqrt(gamma*kb/mu/mp.*tm.*get_tvir)./1e5;  % sonic velocity in km/sec 
% 
% % find shock velocity 
% 
% prat=1.621/0.5;
% u1=0.5.*get_vvir;
% %ushock2 = u1+cs*sqrt(1+(gamma+1)/(2*gamma)*(prat-1))
% mach2= sqrt((prat*(gamma+1)+(gamma-1))/(2*gamma))
% shock2=u1+cs.*mach1
% 
% %% entropy jump
% m=1:0.01:10; %mach number 
% g=gamma;
% gm1=g-1;
% gp1=g+1;
% mm=m.^2;
% sr=(gm1.*mm+2).^g.*(2.*g.*mm-gm1).*gp1.^(-gp1).*mm.^(-g);
