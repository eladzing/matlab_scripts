%% plot to-model for ram pressure stripping

eta=-5:-2; % virial mass ratio between halo and sub


rp=0.01:0.01:3; %position in Halo in units of rv of halo

alfa=0.5;
rs=[];
for i=1:length(eta)
    
rs(end+1,:)=1/sqrt(alfa)*(10.^eta(i))^(1/3)*rp; % stripping radius inunits of rv of sub
end

figure 
h=[];

h(1)=plot(rp,rs(1,:),'-r','linewidth',2.5,...
    'DisplayName','$\mu=10^{-5}$');
hold on
h(2)=plot(rp,rs(2,:),'-b','linewidth',2.5,...
    'DisplayName','$\mu=10^{-4}$');

h(3)=plot(rp,rs(3,:),'color',[0 0.75 0],'linewidth',2.5,...
    'DisplayName','$\mu=10^{-3}$');

h(4)=plot(rp,rs(4,:),'-m','linewidth',2.5,...
    'DisplayName','$\mu=10^{-2}$');


hl=legend(h);
set(hl,'Fontsize',13,'Interpreter','latex','Location','NorthWest')
ylim([0 0.5])
set(gca,'Fontsize',14,'box','on','Ytick',0:0.05:0.5)
grid 
xlabelmine('$r_p/R_c$')
ylabelmine('$\ell_s/R_{sat}$')



