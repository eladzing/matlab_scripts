lw=1.5;
h=[];
figure

beta={'1/2.5' '1/2' '2/3' '1' '1.5' '2' '2.5' 'rand'};
colo=brewermap(8,'Set1');
for i=1:8
    bCen=qFrac(i).bCen;
    qf=qFrac(i).qfProjcnt./qFrac(i).cntProj;
h(i)=plot(bCen,qf,'-','color',colo(i,:),...
    'linewidth',lw,'Displayname',beta{i});
if i==1; hold on;  end
% h(2)=plot(bCen,qf2./cnt2,'rx-',...
%        'linewidth',lw,'Displayname','projected');
end

h(end+1)=plot(Wetzel12_fig15(2,:),Wetzel12_fig15(1,:),'sk',...
'DisplayName','Wetzel+12','markersize',8);

ylim([0 1])
grid

hl=legend(h);
set(hl,'Interpreter','latex','fontsize',14,'Location','NorthEastOutside');

xlabelmine('$r/R_{\mathrm{200,m}}$',16);
ylabelmine('quenched fraction',16);
set(gca,'fontsize',14)



