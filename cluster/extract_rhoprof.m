function [rp lrp rot lrot rog lrog rodm lrodm mto lmto rv mv clname]=extract_rhoprof(ind)
load mat_files/rho_profs_a1.mat
clname=rho_profs{ind,1};
rax=rho_profs{ind,2};
rhot=rho_profs{ind,3}; 
rhodm=rho_profs{ind,4}; 
rhog=rho_profs{ind,5}; 
mtot=rho_profs{ind,6};

rp=rax(1,:); %rho_profs{ind,2};
lrp=rax(2,:); %rho_profs{ind,3};
rot=rhot(1,:);
lrot=rhot(2,:);
rog=rhog(1,:);
lrog=rhog(2,:);
rodm=rhodm(1,:);
lrodm=rhodm(2,:);
mto=mtot(1,:);
lmto=mtot(1,:);


%rot=rho_profs{ind,4};
%lrot=rho_profs{ind,5};
%rodm=rho_profs{ind,6};
%lrodm=rho_profs{ind,7};

virs=rho_profs{ind,7};
rv=virs(1)
lrv=log10(rv)
mv=virs(2);

rp=log10(rp);
lrp=log10(lrp);
rot=log10(rot);
lrot=log10(lrot);
rog=log10(rog);
lrog=log10(lrog);
rodm=log10(rodm)-0.5;
lrodm=log10(lrodm)-0.5;
mto=log10(mto);
lmto=log10(lmto);

%draw density profiles
yl=[min(lrog) 0.5*(min(lrog)+max(lrot))];
lrv=[lrv lrv];
figure;
plot(lrp,lrot,lrp,lrodm,lrp,lrog)
hold on
plot(lrv,yl,'k:');
hold off

%draw mass profile
figure; 
plot(lrp,lmto)
hold on
plot(lrv,[min(lmto) max(lmto)])
plot([min(lrp) max(lrp)],[log10(mv) log10(mv)]);
hold off






