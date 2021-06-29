%% draw Density profiles (gas and total for all halos 

list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
%mask=[ 0   0   1   0   0   1   0   0  0  1  0  0  1  1  1  0];
mask=true(size(list));
syst='win';

arrange_shockedge


for i=1:length(list);
    %continue % THIS IS HERE SINCE THESE PICS ARE DONE
    if~mask(i)
        continue;
    end
    for j=1:2
        switch j
            case 1
                new_env(sprintf('CL%d',list(i)),'csf','a1');
                edge1=edge_a1(i,1:2);
                edge1=edge1(edge1>0);
                edge2=edge_a1(i,3:4);
                edge2=edge2(edge2>0);
            case 2
                new_env(sprintf('CL%d',list(i)),'csf','a06');
                edge1=edge_a06(i,1:2);
                edge1=edge1(edge1>0);
                edge2=edge_a06(i,3:4);
                edge2=edge2(edge2>0);
        end
        
        global DEFUALT_PRINTOUT_DIR;
        global zred  
        global CLUSTER
        global aexpn
        
        r0=find_inner_reslim();
        rv=get_rvir();rv200=get_rvir(200);
        rr=r0:r0:15;rr=rr./(1+zred);
                
        [RHOG_Profile RHOTOT_Profile] = read_RHO_Profiles(rr);
        
        loglog(rr./rv,RHOG_Profile,rr./rv,RHOTOT_Profile,'linewidth',2)
        hold on
        loglog([1 1].*rv200./rv,ylim,'k-.','linewidth',2) 
        for k=1:length(edge1)
            loglog([1 1].*edge1(k)./rv,ylim,'k--','linewidth',2) 
        end
        for k=1:length(edge2)
            loglog([1 1].*edge2(k)./rv,ylim,'k--','linewidth',2) 
        end
        loglog(xlim,[1 1].*rho_mean(zred),'r-.','linewidth',2)
        
        hold off
        grid
        xlim([r0./rv 10])
        
        
        hl=legend('Gas','Total');
        set(hl,'interpreter','latex','Location','SouthWest')
        
        xlabelmine('$r/R_{\mathrm{vir}}$');
        ylabelmine('$\rho\,[\mathrm{M_\odot/Mpc^3}]$');
        
        titlemine(sprintf('%s Density Profiles, $z=%3.2g$',CLUSTER,zred));       
        
        name=sprintf('%s/%s_rho_prof_%s.%s',DEFUALT_PRINTOUT_DIR,CLUSTER,aexpn,'%s');
        %exportfig(gcf,sprintf(name,'png'),'format','png');
        %exportfig(gcf,sprintf(name,'eps'));
    end
end
            
        

