x=0.1:0.1:3.0;
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

global DEFAULT_PRINTUT_DIR

outpath=sprintf('%s/vcm',DEFAULT_PRINTUT_DIR);

v=zeros(length(x),3);
vv=x;
h=0.7;
b=[1 2 4 8];
b=0.5.*b./h;

for i=1:16
    
    cl=sprintf('CL%d',list(i))
    
    new_env(cl,'csf','a06');
       
    vvir=get_vvir;
    rv=get_rvir;
    
    for j=1:20
        v(j,:)=vcmlist(i,:,j)./vvir;
        vv(j)=sqrt(sum(v(j,:).^2));
    end
    %figure
    plot(x,v(:,1),x,v(:,2),x,v(:,3),x,vv,'linewidth',2.0);
    legend('v_x','v_y','v_z','v');legend('boxoff');
    hold on 
    yl=get(gca,'Ylim');
    global aexp
    for k=1:4
        xx=aexp*b(k)./rv;
        if(xx <= 2) 
            xl=[xx xx];
            plot(xl,yl,'--r');
        end
    end
    hold off
    
    xlabelmine('$r/R_{vir}$')
    ylabelmine('$V_{cm}/V_{vir}$')
    set(gca,'Xtick',0:0.1:2);
    
    grid
    titlemine(sprintf('%s $V_{cm}(<r)$, $R_v=%d V_v=%d$',cl,rv,vvir));
    
    exportfig(gcf,sprintf('%s/CL%d_vcm_prof.png',outpath,list(i)),'format','png');
end
    
