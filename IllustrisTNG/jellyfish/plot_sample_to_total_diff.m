
zr50=illustris.utils.snap2redshift(outStruct.snap50);
zr100=illustris.utils.snap2redshift(outStruct.snap100);
%%
myFigure

%%
tt=tiledlayout(2,1);
nexttile
h(1)=semilogx(1+zr50,outStruct.cjfSatNum50./outStruct.totalSatNum50,'-o',...
    'Displayname','TNG50');
hold on
h(2)=semilogx(1+zr100,outStruct.cjfSatNum100./outStruct.totalSatNum100,'-s',...
    'Displayname','TNG100');
legend(h,'Interpreter','latex','fontsize',16,'location','northwest');
ylim([0.5 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabelmine('$N_\mathrm{CJF}/N_\mathrm{Total}$ Sats') ;

nexttile
h2(1)=semilogx(1+zr50,outStruct.cjfhostNum50./outStruct.totalHostNum50,'-o',...
    'Displayname','TNG50');
hold on
h2(2)=semilogx(1+zr100,outStruct.cjfhostNum100./outStruct.totalHostNum100,'-s',...
    'Displayname','TNG100');
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabelmine('$N_\mathrm{CJF}/N_\mathrm{Total}$ Hosts');
ylim([0.85 1])

tt.XLabel.String='$\log(1+z)$';
tt.XLabel.FontSize=18;
tt.XLabel.Interpreter='latex';

% tt.YLabel.String='$N_\mathrm{CJF}/N_\mathrm{Total}$';
% tt.YLabel.FontSize=18;
% tt.YLabel.Interpreter='latex';
tt.TileSpacing='compact';
tt.Padding='compact';
%ylabelmine( o. of CJF sats / No. ofs sats total') 
