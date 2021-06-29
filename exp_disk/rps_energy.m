%% units  are M=Msun L=kpc t=sec.

nn=300;
rlim=[-2 1];
etap_lim=[-2 0.5];

rl=rlim(1):diff(rlim)/(nn-1):rlim(2);
el=etap_lim(1):diff(etap_lim)/(nn-1):etap_lim(2);

r=10.^rl;


etap=10.^el;
[rg eg]=meshgrid(r,etap);

%r=0.01:0.01:10;
units;

gamma=0.1;
alpha=1;
fc=0.2; 
Mc=1e16;
%etap=0.1:0.01:2;

fg=0.2;
beta=0.5;
fb=0.05;

Rd=10;

Gn=G.*(kpc^3/Ms);

sigma_s=10^(7);


poteng=expdisk_pot(rg,'fg',fg,'beta',beta)+bulge_potential(rg,'fb',fb,'point');
gacc=expdisk_accel(rg,'fg',fg,'beta',beta)+bulge_accel(rg,'fb',fb,'point');
rhogas=expdisk_density(rg,'gas','fg',fg,'beta',beta);

poteng=poteng.*pi.*sigma_s.*Rd.*Gn;
gacc=gacc.*pi.*Gn.*sigma_s;
rhogas=rhogas.*sigma_s;

tdyn=gamma*2*pi*sqrt(Rd.*rg./gacc);

lhs=poteng.*rhogas./tdyn;


rhs=(alpha*fc*Gn^(3/2)/8/pi*(4*pi*rho_mean(0)/1e9*deltavir(0))^(1/3)).*Mc.*eg.^(-2);

mat=log10(lhs./rhs);
md=exp_disk_mass(r,beta);
%% plot 
ff=3;
figure(ff)
imagesc(rl,el,mat)
set(gca,'Ydir','normal','Fontsize',12);
bar=colorbar;
%set(get(bar,'Title'),'String','$\%M_{strip}$','Fontsize',12,'Interpreter','latex');

xlabelmine('$\log(R_{disk})\,[\mathrm{R_d}]$')
ylabelmine('$\log(\eta_{p})\,[\mathrm{R_c}]$')

hold on
cc=0.000000001:1;
C=contour(rl,el,mat,cc,'-k','linewidth',1.5);
%clabel(C,'color','white')
hold off