%% test effect of point mass approx. on bulge contribution 

eta=0.1:0.01:100;
fg=0.2;
beta=0.5;
fb=[0.1 0.25] ;
xi=[3 5 10 20];
f0=bfunc(eta,1)+fg.*beta^3.*bfunc(eta,beta);

fbp1=2*fb(1)./eta.^3;

fb11=2*fb(1)*xi(1)^2./(eta.*(1+xi(1).*eta).^2);
fb12=2*fb(1)*xi(2)^2./(eta.*(1+xi(2).*eta).^2);
fb13=2*fb(1)*xi(3)^2./(eta.*(1+xi(3).*eta).^2);
fb14=2*fb(1)*xi(4)^2./(eta.*(1+xi(4).*eta).^2);

fbp2=2*fb(2)./eta.^3;

fb21=2*fb(2)*xi(1)^2./(eta.*(1+xi(1).*eta).^2);
fb22=2*fb(2)*xi(2)^2./(eta.*(1+xi(2).*eta).^2);
fb23=2*fb(2)*xi(3)^2./(eta.*(1+xi(3).*eta).^2);
fb24=2*fb(2)*xi(4)^2./(eta.*(1+xi(4).*eta).^2);



figure 
h=[];

h(1)=loglog(eta,abs((f0+fbp1)./(f0+fb11)-1),'-b','Displayname','$\xi=3,f_b=0.1$','linewidth',2);
hold on
h(2)=loglog(eta,abs((f0+fbp1)./(f0+fb12)-1),'-r','Displayname','$\xi=5,f_b=0.1$','linewidth',2);
h(3)=loglog(eta,abs((f0+fbp1)./(f0+fb13)-1),'-g','Displayname','$\xi=10,f_b=0.1$','linewidth',2);
h(4)=loglog(eta,abs((f0+fbp1)./(f0+fb14)-1),'-m','Displayname','$\xi=20,f_b=0.1$','linewidth',2);

h(5)=loglog(eta,abs((f0+fbp2)./(f0+fb21)-1),'--b','Displayname','$\xi=3,f_b=0.25$','linewidth',2);
h(6)=loglog(eta,abs((f0+fbp2)./(f0+fb22)-1),'--r','Displayname','$\xi=5,f_b=0.25$','linewidth',2);
h(7)=loglog(eta,abs((f0+fbp2)./(f0+fb23)-1),'--g','Displayname','$\xi=10,f_b=0.25$','linewidth',2);
h(8)=loglog(eta,abs((f0+fbp2)./(f0+fb24)-1),'--m','Displayname','$\xi=20,f_b=0.25$','linewidth',2);

ylim([1e-3 1e2])
xlim([0.1 50])
hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14);
grid
xlabelmine('$\eta$')
ylabelmine('$\delta F/F$')