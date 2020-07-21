eta=logspace(-3,6,10000);
alfa=[1 1.2 1.4 1.5];

mf=(1-exp(-eta).*(1+eta)).*100;
h=[];
figure
for i=1:length(alfa)
    
    sfr=(1-exp(-alfa(i).*eta).*(1+alfa(i).*eta)).*100;

    h(i)=plot(mf,sfr,'linewidth',2,'DisplayName',sprintf('$\\alpha=%s$',num2str(alfa(i))));
    hold on
end

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthWest')
grid
xlabelmine('Remaining Mass $[\%]$')
ylabelmine('SFR Reduction $[\%]$');
set(gca,'Fontsize',14,'box','on')