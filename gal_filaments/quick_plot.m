

for i=1:3
    for j=1:3
        figure;imagesc(log10(squeeze(sum(ro(i).Grid,j))));set(gca,'Ydir','Normal');
        titlemine(sprintf('%s $\\rho$ $a=%s$ along %i',fixstring(GAL_NAME),num2str(ro(i).aexp,j)))
        printout_fig(gcf,sprintf('%s_rho_%i_a%s',GAL_NAME,j,num2str(ro(i).aexp)),'png','nopdf')
        
        figure;imagesc(log10(squeeze(mean(tt(i).Grid,j))));set(gca,'Ydir','Normal');
        titlemine(sprintf('%s $T$ $a=%s$ along %i',fixstring(GAL_NAME),num2str(tt(i).aexp,j)))
        printout_fig(gcf,sprintf('%s_tmp_%i_a%s',GAL_NAME,j,num2str(tt(i).aexp)),'png','nopdf')
        
        
    end
end