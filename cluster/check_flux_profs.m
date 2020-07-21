list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%list2=[103 105 106 107];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout';
%streamfluxfile='%s/%s_streamflux_z0_b%d.dat';


plotflag='plot';
pflag='noprint';

vzone=[0.02 0.2];
r_dens=0.01;
r_md=0.2;


for id=1:length(list1)
    clustername=sprintf('CL%d',list1(id))
    new_env(clustername);    
    %stack{id,1}=clustername;
    switch list1(id)
        case {101,102,103,104,105,106,107,5}
            smallbox=2;bigbox=8;
        otherwise
            smallbox=1;bigbox=8;
    end

    for ii=log2(smallbox):log2(bigbox)
        boxx=2^ii;
        [m_dot roprof vrprof r_prof t_prof s_prof rvir mvir tvir vvir]= full_profs(boxx,boxx,'vcm',[0.02 0.2]);
        %[m_dot m_in m_hlfvr roprof r_prof rvir mvir vvir]= penen_profs(boxx,boxx,'vcm',vzone);
    
        [Mg200 ms200 M_dm200]=read_Mass_Profiles(rvir);  
        Mtot= read_MTOT_Profile(r_prof);
        fluxnorm=0.056.*(mvir./1e13).^0.15;
        rhovir=3.*mvir./(4.*pi.*rvir.^3);       
        
        %get stream flux
        if(~strcmp(clustername,'CL5'))
               [rps,fls]=get_streamflux(clustername,boxx,boxx);   
            flxs(:,ii)=(fls./Mg200)./fluxnorm.*1e9;
            rpss(:,ii)=rps./rvir;
        end
        
        flux(:,ii)= (m_dot./Mg200)./fluxnorm.*1e9;
        %flxs(:,ii)=(fls./Mg200)./fluxnorm.*1e9;
        rp(:,ii)=r_prof./rvir;
        %rpss(:,ii)=rps./rvir;
        
        rop(:,ii)=roprof./rhovir;
        rhovir=3.*mvir./(4.*pi.*rvir.^3);   
        
        tpr(:,ii)=t_prof./tvir;
        spr(:,ii)=s_prof;
        vrp(:,ii)=vrprof./vvir;
        
    end
        if strcmp(plotflag,'plot')
        figure; 
   
         subplot(3,2,1);
    semilogx(rp,flux);grid;xlim([1e-3 10]);
    %xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('$\dot{\mathrm{M}}$','Fontsize',12,'Interpreter','latex');
    ylabel('total flux','Fontsize',12);
    %title(sprintf('Total Flux, $\dot{mathrm{M}}_{vir}=%s \mathrm{Gyr}^{-1}$ ', num2str((fluxnorm.*1e9),'%1.2d')),'Fontsize',12,'Interpreter','latex');
    title(sprintf('%s - multi box profiles',clustername),'Fontsize',12);
    set(gca,'XTick', [1e-3 1e-2 1e-1 1 2.5],'Fontsize',12);%,'Interpreter','latex');
    
    if(~strcmp(clustername,'CL5'))
     subplot(3,2,2);
    semilogx(rpss,flxs);grid;xlim([1e-3 10]);
    %xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('$\dot{\mathrm{M}}$','Fontsize',12,'Interpreter','latex');
    ylabel('Streamline Flux','Fontsize',12);
    %title(sprintf('Streamline Flux, $\dot{mathrm{M}}_{vir}=%s \mathrm{Gyr}^{-1}$ ', num2str((fluxnorm.*1e9),'%1.2d')),'Fontsize',12,'Interpreter','latex');
    set(gca,'XTick', [1e-3 1e-2 1e-1 1 2.5],'Fontsize',12);%,'Interpreter','latex');   
    end
    subplot(3,2,3);
    loglog(rp,rop);grid;xlim([1e-3 10]);
    %xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('$\rho/\rho_{vir}','Fontsize',12,'Interpreter','latex');
    ylabel('\rho/\rho_{vir}','Fontsize',12);
    %title(sprintf('Density, $\rho_{vir}=%s \mathrm{M}_{\odot}/\mathrm{Mpc}^3$', num2str(rhovir,'%1.2d')),'Fontsize',12,'Interpreter','latex');
    set(gca,'XTick', [1e-3 1e-2 1e-1 1 2.5],'Fontsize',12);%,,'Interpreter','latex');
    
     subplot(3,2,4);
    loglog(rp,tpr);grid;xlim([1e-3 10]);
    %xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('$T/T_{vir}$','Fontsize',12,'Interpreter','latex');
    ylabel('T/T_{vir}','Fontsize',12);
    %title(sprintf('Temperature, T_{vir}=%s K', num2str(tvir,'%1.2d')),'Fontsize',12,'Interpreter','latex');
    set(gca,'XTick', [1e-3 1e-2 1e-1 1 2.5],'Fontsize',12);%,,'Interpreter','latex');
    
    subplot(3,2,5);
    loglog(rp,spr);grid;xlim([1e-3 10]);
    %xlabel('$r/R_{vir}$','Fontsize',12,'Interpreter','latex');
    xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('S','Fontsize',12,'Interpreter','latex');
    ylabel('S','Fontsize',12);
    %title('Entropy','Fontsize',12,'Interpreter','latex');
    set(gca,'XTick', [1e-3 1e-2 1e-1 1 2.5],'Fontsize',12);%,'Interpreter','latex');
    
     subplot(3,2,6);
    semilogx(rp,vrp);grid;xlim([1e-3 10]);
    %xlabel('$r/R_{vir}$','Fontsize',12,'Interpreter','latex');
    xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('$V_r/V_{vir}$','Fontsize',12,'Interpreter','latex');
    ylabel('V_r/V_{vir}','Fontsize',12);
    %title(sprintf('Radial Velocity, $V_{vir}=%s \mathrm{km/sec}', num2str(vvir,'%1.2d')),'Fontsize',12,'Interpreter','latex');
    set(gca,'XTick', [1e-3 1e-2 1e-1 1 2.5],'Fontsize',12);%,'Interpreter','latex');
    
    if strcmp(pflag,'print')
        saveas(gcf,sprintf('%s/%s_all_profs.png',result_dir,clustername));
    end
   
        end
end

        