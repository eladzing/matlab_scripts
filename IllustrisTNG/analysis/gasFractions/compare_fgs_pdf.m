

global DEFAULT_MATFILE_DIR


load([DEFAULT_MATFILE_DIR '/fgs_stellarMass_gaussianFit_byBin_snp_99_TNG300.mat'])



binEdge=[9 9.2500 9.5000 9.7500 10 10.2500 10.5000 10.7500 11 11.2500 11.5000];
be=9:0.25:11.5;


fgs=-5:0.001:1.5;
hf=figure('color','w','position',[83 433 1710 883]);




for i=1:length(ft(1).mu1)
    h=[];
    
    pdf1=zeros(length(ft),length(fgs));
    for k=1:length(ft)
        pdf=ft(k).a1(i).*exp(-0.5.*((fgs-ft(k).mu1(i))./ft(k).sig1(i)).^2) + ...
            ft(k).a2(i).*exp(-0.5.*((fgs-ft(k).mu2(i))./ft(k).sig2(i)).^2);
        pdf1(k,:)=pdf./trapz(fgs,pdf);
        
    end
    
    %figure(hf)
    subplot(2,5,i)
    ax=gca;
    
    
    for k=1:length(ft)
        tag=['hist bins=' num2str(ft(k).histBinNum)];
        h(k)=plot(fgs,pdf1(k,:),'Displayname',tag);
        hold on
    end
    
    if i==1
        hl=legend(h);
        set(hl,'Interpreter','latex','fontsize',14,'Location','NorthWest');
    end
    
    %ax.Legend.Visible='off';
    
    grid
    xlim([-5 1.5])
    if i>5
        xlabelmine( '$\log\,f_\mathrm{gs}$');
    else
        ax.XLabel.Visible='off';
        
    end
    if i==1 || i==6
        ylabelmine( '$PDF$');
    else
        ax.YLabel.Visible='off';
        
    end
    titlemine([ num2str(be(i)) '-' num2str(be(i+1)) ]);
    
end
    