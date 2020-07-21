
outpath='~/work/clusters/printout';
cflist=[104 105 107 3 6 7 10 14];
ncflist=[101 102 103 106 5 9 11 24];

r_flist=[104 3 5 7 10 14]; %relaxed
ur_flist=[101 102 103 105 106 107 6 9 11 24];

add=5;
for i=1:14
    name=thermal_eng_stack{i,1}
    rb=thermal_eng_stack{i,5};

    eth=thermal_eng_stack{i,6};
    pdv=thermal_eng_stack{i,7};

    pars=thermal_eng_stack{i,8};
    rv=pars(2);
    res=coolprofs{i,6};
    qb=res(2,:); clear res;
    res=coolprofs{i,7};
    ind=res(1)+add; clear res;
        
    cq1=fliplr(cumsum(fliplr(qb(ind:end))));
    cp1=fliplr(cumsum(fliplr(pdv(ind:end))));
    cet1=fliplr(cumsum(fliplr(eth(ind:end))));
    rp1=rb(ind:end)./rv;
    rat1=(cp1-cq1)./cet1.*1e9;
    clear ind
    
    ind=find(rb>0.05*rv,1,'first');
    cq2=fliplr(cumsum(fliplr(qb(ind:end))));
    cp2=fliplr(cumsum(fliplr(pdv(ind:end))));
    cet2=fliplr(cumsum(fliplr(eth(ind:end))));
    rp2=rb(ind:end)./rv;
    rat2=(cp2-cq2)./cet2.*1e9;
    
    
    for j=1:length(cflist)
        if(strcmp(name,sprintf('CL%d',cflist(j))))
            tag1='CF'
            break
        end
    end

    for j=1:length(ncflist)
        if(strcmp(name,sprintf('CL%d',ncflist(j))))
            tag1='NCF'
            break
        end
    end
    for j=1:length(r_flist)
        if(strcmp(name,sprintf('CL%d',r_flist(j))))
            tag2='RLX'
            break
        end
    end
    for j=1:length(ur_flist)
        if(strcmp(name,sprintf('CL%d',ur_flist(j))))
            tag2='URLX'
            break
        end
    end
    
    figure
    
    semilogx(rp1,rat1,rp2,rat2,'Linewidth',2)
    if(strcmp(name,'CL101'))
        legend('r_{min}=inner','r_{min}=0.05R_{vir}','Location','SouthEast')
    else
        legend('r_{min}=inner','r_{min}=0.05R_{vir}')
    end
   
    xlabel('r/R_{vir}','Fontsize',14)
    ylabel('$$\dot{E}_{th}/E_{th} [1/Gyr] $$','Interpreter','latex','Fontsize',14)
    set(gca,'Fontsize',14)
    title(sprintf('%s (%s,%s) E_{th} rate',name,tag2,tag1),'fontsize',14)
    grid
    saveas(gcf,sprintf('%s/%s_ethrate_flip.png',outpath,name));

    
    figure
    semilogx(rp1,cp1,rp2,cp2,rp1,cq1,rp2,cq2,'Linewidth',2)
    legend('pdv_1','pdv_2','q_1','q_1')
    xlabel('r/R_{vir}','Fontsize',14)
    ylabel('$$\dot{E}_{th}  $$','Interpreter','latex','Fontsize',14)
    title(sprintf(' %s (%s,%s) - sink/source',name,tag2,tag1),'fontsize',14)
    grid
    saveas(gcf,sprintf('%s/%s_ethdot_flip.png',outpath,name));
    
    
    figure
    semilogx(rp1,cet1,rp2,cet2,'Linewidth',2)
    legend('E_{th} - 1','E_{th} - 2')
    xlabel('r/R_{vir}','Fontsize',14)
    ylabel('$${E}_{th}  $$','Interpreter','latex','Fontsize',14)
    title(sprintf('%s (%s,%s) - E_{th}',name,tag2,tag1),'fontsize',14)
    grid
    saveas(gcf,sprintf('%s/%s_eth_flip.png',outpath,name));
    
   
    
    pause
    
end
close all