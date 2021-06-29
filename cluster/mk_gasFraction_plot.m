%% calculate gas fraction profiles in clusters in rlx / urlx populations
penetration_depth

r0=0.05:0.01:4;
r=r0(1:end-1)+0.5*diff(r0);
rr=cat(2,r,fliplr(r));


fg=zeros(length(clust),length(r));
fb=fg;
h=[];
hl=[];
for k=1:2
    switch k
        case 1
            a='a1';
        case 2
            a='a06';
    end
    
    
    for i=1:length(clust)
        
        new_env(clust(i),a)
        
        [mgas,mstar,mdm] = read_Mass_Profiles(r0.*get_rvir);
        
        fg(i,:)=diff(mgas)./diff(mgas+mstar+mdm);
        fb(i,:)=diff(mgas+mstar)./diff(mgas+mstar+mdm);
    end
    
    
    
    cc=zeros(length(clust),3);
    cc(urlx,1)=1;
    cc(rlx,3)=1;
    
    %% gas fraction plots
    figure
    for i=1:length(clust)
        semilogx(r,fg(i,:),'color',cc(i,:),'linestyle',':')
        hold on
    end
    global zred 
    
    mu=mean(fg(urlx,:));
    stu=std(fg(urlx,:));
    pu=cat(2,mu+0.5.*stu,fliplr(mu-0.5.*stu));
    h(1)=semilogx(r,mu,'r-','linewidth',2,'DisplayName','NR');
    patch(rr,pu,'r','FaceAlpha',0.2)
    
    mr=mean(fg(rlx,:));
    str=std(fg(rlx,:));
    pr=cat(2,mr+0.5.*str,fliplr(mr-0.5.*str));
    h(2)=semilogx(r,mr,'b-','linewidth',2,'DisplayName','R');
    patch(rr,pr,'b','FaceAlpha',0.2)
    
    xlim([4e-2 5])
    
    xlabelmine('$r/R_{\mathrm{vir}}$')
    ylabelmine('$f_{\mathrm{gas}}$')
    grid
    set(gca,'Fontsize',14);

    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
    
    titlemine(sprintf('Gas Fraction $z=%s$',num2str(zred)));
    
    %% baryon fraction plots (gas+stars)
    figure
    for i=1:length(clust)
        semilogx(r,fb(i,:),'color',cc(i,:),'linestyle',':')
        hold on
    end
    
    
    mu=mean(fb(urlx,:));
    stu=std(fb(urlx,:));
    pu=cat(2,mu+0.5.*stu,fliplr(mu-0.5.*stu));
    h(1)=semilogx(r,mu,'r-','linewidth',2,'DisplayName','NR');
    patch(rr,pu,'r','FaceAlpha',0.2)
    
    mr=mean(fb(rlx,:));
    str=std(fb(rlx,:));
    pr=cat(2,mr+0.5.*str,fliplr(mr-0.5.*str));
    h(2)=semilogx(r,mr,'b-','linewidth',2,'DisplayName','R');
    patch(rr,pr,'b','FaceAlpha',0.2)
    
    xlim([4e-2 5])
    
    xlabelmine('$r/R_{\mathrm{vir}}$')
    ylabelmine('$f_{\mathrm{baryon}}$')
    
    grid
    set(gca,'Fontsize',14);
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',14)
     titlemine(sprintf('Baryon Fraction $z=%s$',num2str(zred)));
end