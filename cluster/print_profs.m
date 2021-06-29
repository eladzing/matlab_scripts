list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
typ='csf';% 'adiabatic');
aexp='a1'; %'a06'];
%list2=[103 105 106 107];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout';
%result_dir='/home/eladzing/work/sshfs/titan3/cold_flows/tit3/printout/';
%result_dir='/home/titan3/eladzing/cold_flows/printout';

plotflag='plot';
pflag='noprint';
%hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    
hubblefl='hub';



for id=1:length(list1)
    %clustername=sprintf('CL%d',list1(id))
    %new_env(clustername,typ,aexp);    
    

    %[m_dot roprof vrprof r_prof t_prof s_prof rvir mvir tvir vvir]= full_profs(smallbox,bigbox,'vcm',[0 1],hubblefl);
    %[Mg_RV ms_RV M_dm_RV]=read_Mass_Profiles(rvir);  
    %Mtot= read_MTOT_Profile(r_prof);
    %rhovir=3.*mvir./(4.*pi.*rvir.^3); 
     
    load(sprintf('mat_files/stk_profs_%s.mat',aexp),'stack');
      
    clustername=stack{id,1}  
    rp=stack{id,2}; %=r_prof/rvir; %r
    md=stack{id,3}; %=m_dot./Mg_RV; %Flux profile
    tp=stack{id,4};% =t_prof./tvir; %temperature profile 
    rop=stack{id,5};% =roprof./rhovir; %density profile 
    vrp=stack{id,6};% =vrprof./vvir; %radial velocity profile 
    sp=stack{id,7};% =s_prof./(tvir./rhovir.^(2/3)); %full mass flux
    norms=stack{id,8};% =[mvir,rvir,Mg_RV,tvir,rhovir,vvir,(tvir./rhovir.^(2/3))];
    
    mvir=norms(1);
    rvir=norms(2);
    Mg_RV=norms(3);
    tvir=norms(4);
    rhovir=norms(5);
    vvir=norms(6);
    svir= norms(7);
    
    
    fluxnorm=0.056.*(mvir./1e13).^0.15;
    flux= (md)./fluxnorm.*1e9;
    
    if strcmp(plotflag,'plot')
        figure; 
        
        %subplot(3,2,1);
        % title(sprintf('%s' ,clustername),'Fontsize',12);
        % set(gca, 'XTick', [], 'YTick', []);
        % text(0.1,0.9,sprintf('M_{vir}=%s',num2str(mvir,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
        % text(0.1,0.7,sprintf('M_{gas}=%s',num2str(Mg_RV,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
    
        subplot(3,1,1:2);
        semilogx(rp,flux);grid;xlim([1e-3 3]);ylim([-6 2]);
        %xlabel('r/R_{vir}','Fontsize',12);
        ylabel('dM/dt [virial dm/dt]','Fontsize',12);
        title(sprintf('%s Flux, Virial dM/dt=%s Gyr^{-1},%s',clustername,num2str(fluxnorm,'%1.2d'),aexp),'Fontsize',12);
        line([1e-3 1e1],[0 0],'Color','k');
        set(gca,'Fontsize',12);
    
        subplot(3,1,3)
        semilogx(rp,vrp);grid;xlim([1e-3 3]);ylim([-1 0.25]);
        xlabel('r/R_{vir}','Fontsize',12);
        ylabel(sprintf('V_r/V_{vir} [V_{vir}=%s km/sec]',num2str(vvir,4)),'Fontsize',12);
        %title(sprintf('Radial Velocity, V_{vir}=%s km/sec', num2str(vvir,'%1.2d')),'Fontsize',12);
        set(gca,'Fontsize',12);  
    
        if strcmp(pflag,'print')
            saveas(gcf,sprintf('%s/%s_flux_vr_prof_%s.png',result_dir,clustername,aexp));
        end
    
    
        figure;
        subplot(2,2,1);
        title(sprintf('%s' ,clustername),'Fontsize',12);
        set(gca, 'XTick', [], 'YTick', []);
        text(0.1,0.9,sprintf('M_{vir}=%s M_{sun}',num2str(mvir,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
        text(0.1,0.7,sprintf('M_{gas}=%s M_{sun}',num2str(Mg_RV,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
        text(0.1,0.5,sprintf('R_{vir}=%s Mpc',num2str(rvir,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
        text(0.1,0.3,sprintf('%s',aexp), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
    
    
        subplot(2,2,2);
        loglog(rp,tp);grid;xlim([1e-3 3]);ylim([5e-2 13]);
        %xlabel('r/R_{vir}','Fontsize',12);
        ylabel('T/T_{vir}','Fontsize',12);
        title(sprintf('Temperature, T_{vir}=%s K', num2str(tvir,'%1.2d')),'Fontsize',12);
        set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick', [1e-3 1e-2 1e-1 1 10 100 1000 10000],'Fontsize',12);
    
        subplot(2,2,3);
        loglog(rp,rop);grid;xlim([1e-3 3]);ylim([8e-3 1e4]);
        xlabel('r/R_{vir}','Fontsize',12);
        ylabel('\rho_{gas}/\rho_{vir}','Fontsize',12);
        title(sprintf('Gas Density, %s{vir}=%s M_{sun}/Mpc^3', texlabel('rho'),num2str(rhovir,'%1.2d')),'Fontsize',12);
        set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick', [1e-3 1e-2 1e-1 1 10 100 1000 10000],'Fontsize',12);
    
    %subplot(2,2,4);
    %semilogx(rp,vrprof./vvir);grid;xlim([1e-3 3]);%ylim([ ]);
    %xlabel('r/R_{vir}','Fontsize',12);
    %ylabel('V_r/V_{vir}','Fontsize',12);
    %title(sprintf('Radial Velocity, V_{vir}=%s km/sec', num2str(vvir,'%1.2d')),'Fontsize',12);
    %set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick', [1e-3 1e-2 1e-1 1 10],'Fontsize',12);  
    
        subplot(2,2,4);
        loglog(rp,sp);grid;xlim([1e-3 3]);ylim([1e-4 1e2]);
        xlabel('r/R_{vir}','Fontsize',12);
        ylabel('S/S_{vir}','Fontsize',12);
        title(sprintf('Entropy, S_{vir}=%s K Mpc^2 /{M_{sun}}^{2/3}',num2str(svir,'%1.2d')),'Fontsize',12);
        set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick', [1e-3 1e-2 1e-1 1 10 100 1000 10000],'Fontsize',12);
    
        if strcmp(pflag,'print')
            saveas(gcf,sprintf('%s/%s_tros_profs_%s.png',result_dir,clustername,aexp));
        end
    
    end
        
    %close all
end


  
