cl=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
vf=[1.39 1.681 1.239 1.527 2.172 2.009 2.624 2.87 3.037 2.363 1.95 2.954];

cc= brewermap(8,'Set1');
for i=1:1%length(cl)
    
    new_env(cl(i))
    global hub;
    rv=get_rvir();
    h=[];
    il=1;
    figure
    for boxx=[1 2 4 8]
        
        [hubX, hubY, hubZ] = hubble_flow(boxx);
        Vxx = Vx(boxx)+hubX;
        Vyy = Vy(boxx)+hubY;
        Vzz = Vz(boxx)+hubZ;
        ro=RHOG(boxx);
        rop=MAKE_PROFILE_FROM_CUBE(ro);
        rop2=cumsum(rop);
        vx=MAKE_PROFILE_FROM_CUBE(Vxx.*ro);
        vy=MAKE_PROFILE_FROM_CUBE(Vyy.*ro);
        vz=MAKE_PROFILE_FROM_CUBE(Vzz.*ro);
        vx2=cumsum(vx);
        vy2=cumsum(vy);
        vz2=cumsum(vz);
        vxc=vx2./rop2;
        vyc=vy2./rop2;
        vzc=vz2./rop2;
        vc=sqrt(vxc.^2+vyc.^2+vzc.^2);
        rr=linspace(boxx./2./256,boxx./2,256);
        rr=rr./hub./rv;
        ll=il:256;
        
        h(1)=plot(rr(ll),vxc(ll),'color',cc(1,:),'linewidth',1.8,'DisplayName','$V_x$');%sprintf('$V_x,b%s$',num2str(boxx)));
        hold on
        h(2)=plot(rr(ll),vyc(ll),'color',cc(2,:),'linewidth',1.8,'DisplayName','$V_y$');%sprintf('$V_y,b%s$',num2str(boxx)));
        h(3)=plot(rr(ll),vzc(ll),'color',cc(3,:),'linewidth',1.8,'DisplayName','$V_z$');%sprintf('$V_z,b%s$',num2str(boxx)));
        h(4)=plot(rr(ll),vc(ll),'color',cc(4,:),'linewidth',1.8,'DisplayName','$|V|$');%sprintf('$V_c,b%s$',num2str(boxx)));
        
        if boxx==1
            il=129;
        end
    end
    
    %hl=clickableLegend(h);
    plot(vf(i).*[1 1],ylim,'--k','linewidth',2) 
    hl=legend(h);
    set(gca,'box','on','Fontsize',14)
    set(hl,'Interpreter','latex','Fontsize',14)
    xlabelmine('$r/R_{\mathrm{vir}}$')
    ylabelmine('$V\,[\mathrm{km\,sec^{-1}}]$')
    titlemine(sprintf('CL%s',num2str(cl(i))));
    fname=sprintf('cl%s_vcmProf',num2str(cl(i)));
    printout_fig(gcf,fname);
end