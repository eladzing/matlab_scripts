for i=1:length(ind)
    
    if centralHist(ind(i)).stellarMass(1)<1e10
        
        ssfr=centralHist(ind(i)).sfr./centralHist(ind(i)).stellarMass;
        rad=centralHist(ind(i)).radiusToHost./centralHist(ind(i)).hostR200;
        mm=rad>0.01;
        
        yyaxis left
        plot(centralHist(ind(i)).zred,centralHist(ind(i)).isCentral,'b-',...
            centralHist(ind(i)).zred(mm),rad(mm),'mx')
        ylim([-0.02 2])
        
        yyaxis right
        semilogy(centralHist(ind(i)).zred,ssfr,'-')
        ylabelmine('sSFR');
        
        xlim([-0.02 4])
        xlabel('z')
        titlemine(sprintf('%i %f',ind(i)-1,log10(centralHist(ind(i)).stellarMass(1))));
        
        pause
    end
end