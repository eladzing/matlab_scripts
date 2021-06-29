function res=effective_rad(eff)

%% function for finding the effective radius for a given mass
% fraction of an exponential disk in units of rd

%x=0:0.000001:10;
x=0:1e-4:10;
y=1-exp(-x).*(1+x);

res=interp1(y,x,eff);
%ind=find(y>eff,1,'first');
%res(2)=x(ind);

end




