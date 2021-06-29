
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];



v=zeros(size(vcmlist2,3),3);
vv=1:size(vcmlist2,3);
h=0.7;
b=[1 2 4 8];
b=0.5.*b./h;

for i=1:length(list)
    
    cl=sprintf('CL%d',list(i))
    
    new_env(cl,'csf','a06');
       
    global DEFAULT_PRINTOUT_DIR
    global aexp
    global zred
    global aexpn
    outpath=sprintf('%s/vcm',DEFAULT_PRINTOUT_DIR);
    
    vvir=get_vvir;
    rv=get_rvir;
    rr=squeeze(vcmlist2(i,1,:));
    for j=1:size(vcmlist2,3)
        v(j,:)=vcmlist2(i,2:4,j);%./vvir;
        vv(j)=sqrt(sum(v(j,:).^2));
    end
    figure
    plot(rr,v(:,1),rr,v(:,2),rr,v(:,3),rr,vv,'linewidth',2.0);
    hl=legend('$v_x$','$v_y$','$v_z$','$v$');legend('boxoff');
    set(hl,'Interpreter','latex');
     hold on 
     yl=get(gca,'Ylim');
     for k=1:4
         xx=aexp*b(k)./rv;
         %if(xx <= 2) 
             xl=[xx xx];
             plot(xl,yl,'--r');
         %end
     end
     hold off
    
    xlabelmine('$r/R_{vir}$');
    ylabelmine('$V_{cm}$');
    %set(gca,'Xtick',0:0.1:2);
    
    grid
    
    titlemine(sprintf('%s $V_{cm}(<r)$  $z=%3.2g$ , $R_{vir}=%3.3f $ $ V_{vir}=%3.1f$',cl,zred,rv,vvir));
    exportfig(gcf,sprintf('%s/CL%d_vcm_prof_%s.png',outpath,list(i),aexpn),'format','png');
    %saveas(gcf,sprintf('%s/CL%d_vcm_prof2.png',outpath,list(i)));
    %pause
    
end

