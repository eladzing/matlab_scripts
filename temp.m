figure
hN=[];
hI=[];
for j=1:length(ms)
    dnTag=sprintf('$m_{\\mathrm{sat}}=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(floor(log10(ms(j)))));
    [ff1, xx1]=derive1(msN(:,j),rp');
    [ff2, xx2]=derive1(msI(:,j),rp');
    %hN(j)=plot(rp,msN(:,j),'-','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    hN(j)=plot(xx1, ff1,'-','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    hold on
    hI(j)=plot(xx2,ff2,'--','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    
end
    
grid
xlim([0.1 3])
ylim([1e-2 0.8])

%hl=legend(cat(2,hN,hI));
hl=legend(hN);

set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);

set(gca,'Fontsize',12);
xlabelmine('$r_p/R_c$')
ylabelmine('$\mathrm{d} m/ \mathrm{d} r$')
