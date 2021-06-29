
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
outpath='~/work/clusters/printout';

v=zeros(size(vcmlist2,3),3);
vv=1:size(vcmlist2,3);
h=0.7;
b=[1 2 4 8];
b=0.5.*b./h;

for i=1:16
    
    cl=sprintf('CL%d',list(i))
    
    new_env(cl,'csf','a1');
       
    vvir=get_vvir;
    rv=get_rvir;
    rr=squeeze(vcmlist2(i,1,:));
    for j=1:size(vcmlist2,3)
        v(j,:)=vcmlist2(i,2:4,j);%./vvir;
        vv(j)=sqrt(sum(v(j,:).^2));
    end
    %figure
    dr=diff(rr(1:end-1))+diff(rr(2:end));
    dv=diff(v(1:end-1,:))+diff(v(2:end,:));
    dvv=(diff(vv(1:end-1))+diff(vv(2:end)));
    
    dvv=dvv./abs(vv(2:end-1));
    for j=1:3
        dv(:,j)=dv(:,j)./abs(v(2:end-1,j));
    end
    
    rdr=rr(2:end-1);
    plot(rdr,dv(:,1)./dr,rdr,dv(:,2)./dr,rdr,dv(:,3)./dr,rdr,dvv'./dr,'linewidth',2.0);
    legend('v_x','v_y','v_z','v');legend('boxoff');
     hold on 
     yl=get(gca,'Ylim');
     for k=1:4
         xx=b(k)./rv;
         %if(xx <= 2) 
             xl=[xx xx];
             plot(xl,yl,'--r');
         %end
     end
     hold off
    
    xlabel('$r/R_{vir}$','Fontsize',12,'Interpreter','latex');        
    ylabel('$dv/dr$','Fontsize',12,'Interpreter','latex');
    %set(gca,'Xtick',0:0.1:2);
    
    grid
    title(sprintf('%s dV_{cm}(<r), R_v=%d V_v=%d',cl,rv,vvir));
    
    saveas(gcf,sprintf('%s/CL%d_vcm_prof4.png',outpath,list(i)));
    %pause
    
end

