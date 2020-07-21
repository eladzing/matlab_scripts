%% centrals 
illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBird./countA,'tit',sprintf('All: %s',num2str(countA)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_All_tng100')


illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdC./countC,'tit',sprintf('Centrals: %s',num2str(countC)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Centrals_tng100')


illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdQC./countQC,'tit',sprintf('Quenched Centrals: %s',num2str(countQC)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Centrals_Quenched_tng100')

illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdTC./countTC,'tit',sprintf('trans Centrals: %s',num2str(countTC)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Centrals_trans_tng100')


illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdSFC./countSFC,'tit',sprintf('SF Centrals: %s',num2str(countSFC)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Centrals_SF_tng100')


%% sats
illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdS./countS,'tit',sprintf('Sats: %s',num2str(countS)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Sats_tng100')


illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdQS./countQS,'tit',sprintf('Quenched Sats: %s',num2str(countQS)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Sats_Quenched_tng100')

illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdTS./countTS,'tit',sprintf('trans Sats: %s',num2str(countTS)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Sats_trans_tng100')


illustris.plots.plot_phaseDiagram(xxlim,yylim,cgmBirdSFS./countSFS,'tit',sprintf('SF Sats: %s',num2str(countSFS)),'lines','clims',[-7 0])
printout_fig(gcf,'cgmBird_Sats_SF_tng100')


