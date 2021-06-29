list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
typ='csf';% 'adiabatic');
aexp='a1'; 

%list1=[101 102 103 104 105 106 107 3 6 7 9 10 11 14 24];

%list2=[103 105 106 107];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout';
%streamfluxfile='%s/%s_streamflux_z0_b%d.dat';

%result_dir='/home/eladzing/work/sshfs/titan3/cold_flows/tit3/printout/';
%result_dir='/home/titan3/eladzing/cold_flows/printout';

plotflag='noplot';
pflag='noprint';
hubblefl='hub';

%hpath='/home/alf/eladzing/data/kravtsov/clusters/CL%d';    
vzone=[0.02 0.2];
r_dens=0.01;
r_md=0.2;

% ropar=[];
% flxpar1=[];
% flxpar2=[];
% flxpar3=[];
% flxpar4=[];
% cname=[];

for id=1:length(list1)
    clustername=sprintf('CL%d',list1(id))
    new_env(clustername,typ,aexp);       
    %stack{id,1}=clustername;
%     switch list1(id)
%         case {101,102,103,104,105,106,107,5}
%             smallbox=2;bigbox=2;
%         otherwise
%             smallbox=1;bigbox=2;
%     end

    %[m_dot roprof vrprof r_prof t_prof s_prof rvir mvir tvir vvir]= full_profs(smallbox,bigbox,'vcm',[0.02 0.2]);
    [m_dot m_in m_hlfvr roprof vrprof vrin r_prof rvir mvir vvir]= penen_profs(smallbox,bigbox,'vcm',vzone,hubblefl);
    
    [Mg_RV ms_RV M_dm_RV]=read_Mass_Profiles(rvir);  
    Mtot= read_MTOT_Profile(r_prof);
    rhovir=3.*mvir./(4.*pi.*rvir.^3); 
    
    %get stream flux
%     if(~strcmp(clustername,'CL5'))
%         [rps,fls]=get_streamflux(clustername,smallbox,bigbox);
%     else
%         rps=zeros(size(fls));
%         fls=zeros(size(rps));
%     end
    
    %stack{id,2}=r_prof/rvir; %r
    %stack{id,3}=m_dot./Mg_RV; %Flux profile
    %stack{id,4}=t_prof./tvir; %temperature profile 
    %stack{id,5}=roprof./rhovir; %density profile 
    %stack{id,6}=vrprof./vvir; %radial velocity profile 
    %stack{id,7}=s_prof./(tvir./rhovir.^(2/3)); %full mass flux
    %stack{id,8}=[mvir,rvir,Mg_RV,tvir,rhovir,vvir,(tvir./rhovir.^(2/3))];
    
    
    %fluxnorm=0.056.*(mvir./1e13).^0.15;
    flux1= (m_dot./Mg_RV);%./fluxnorm.*1e9;
    flux2= (m_in./Mg_RV);%./fluxnorm.*1e9;
    flux3= (m_hlfvr./Mg_RV);%./fluxnorm.*1e9;
    %flux4= (fls./Mg_RV);%./fluxnorm.*1e9;
    rp=r_prof./rvir;
    rps=rps./rvir;
    rop=roprof./rhovir;
    
    stackpen{id,1}=clustername;
    stackpen{id,2}=rp;
    stackpen{id,3}=flux1;
    stackpen{id,4}=flux2;
    stackpen{id,5}=flux3;
    stackpen{id,6}=0 %flux4;
    stackpen{id,7}=rps;
    stackpen{id,8}=rop;
    stackpen{id,9}=[mvir,rvir,Mg_RV,rhovir,vvir];
    stackpen{id,10}=vrprof./vvir;
    stackpen{id,11}=vrin./vvir;
end
    save(sprintf('mat_files/stk_penen_%s.mat',aexp),'stackpen');

    
    
    
    
%     ropar(end+1)=interp1(rp,rop,r_dens,'linear');
%     flxpar1(end+1)=interp1(rp,flux1,r_md,'linear');
%     flxpar2(end+1)=interp1(rp,flux2,r_md,'linear');
%     flxpar3(end+1)=interp1(rp,flux3,r_md,'linear');
%     if(~strcmp(clustername,'CL5'))
%         flxpar4(end+1)=interp1(rps,flux4,r_md,'linear');
%     else
%         flxpar4(end+1)=0;
%     end
%     
%     cname{end+1}=clustername;
%     
%     
%     
% end
%     if strcmp(plotflag,'plot')
%     
%         figure;plotmatrix(flxpar',log10(ropar)');
        
    
    
    %subplot(3,2,1);
    % title(sprintf('%s' ,clustername),'Fontsize',12);
    % set(gca, 'XTick', [], 'YTick', []);
    % text(0.1,0.9,sprintf('M_{vir}=%s',num2str(mvir,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
    % text(0.1,0.7,sprintf('M_{gas}=%s',num2str(Mg_RV,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
%     
%     subplot(3,1,1:2);
%     semilogx(r_prof/rvir,flux);grid;xlim([1e-3 3]);ylim([-6 2]);
%     %xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('dM/dt [virial dm/dt]','Fontsize',12);
%     title(sprintf('%s Flux, Virial dM/dt=%s Gyr^{-1}',clustername,num2str(fluxnorm,'%1.2d')),'Fontsize',12);
%     line([1e-3 1e1],[0 0],'Color','k');
%     set(gca,'Fontsize',12);
%     
%     subplot(3,1,3)
%     semilogx(r_prof/rvir,vrprof./vvir);grid;xlim([1e-3 3]);ylim([-1 0.25]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel(sprintf('V_r/V_{vir} [V_{vir}=%s km/sec]',num2str(vvir,4)),'Fontsize',12);
%     %title(sprintf('Radial Velocity, V_{vir}=%s km/sec', num2str(vvir,'%1.2d')),'Fontsize',12);
%     set(gca,'Fontsize',12);  
%     
%     if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_vrflux_prof.png',result_dir,clustername));
%     end
%     
%     
%     figure;
%     subplot(2,2,1);
%      title(sprintf('%s' ,clustername),'Fontsize',12);
%      set(gca, 'XTick', [], 'YTick', []);
%      text(0.1,0.9,sprintf('M_{vir}=%s M_{sun}',num2str(mvir,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
%      text(0.1,0.7,sprintf('M_{gas}=%s M_{sun}',num2str(Mg_RV,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
%      text(0.1,0.5,sprintf('R_{vir}=%s Mpc',num2str(rvir,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
%      
%     
%     
%     subplot(2,2,2);
%     loglog(r_prof/rvir,t_prof./tvir );grid;xlim([1e-3 3]);ylim([5e-2 6]);
%     %xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('T/T_{vir}','Fontsize',12);
%     title(sprintf('Temperature, T_{vir}=%s K', num2str(tvir,'%1.2d')),'Fontsize',12);
%     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);
%     
%     subplot(2,2,3);
%     loglog(r_prof/rvir,roprof./rhovir );grid;xlim([1e-3 3]);ylim([8e-3 1e4]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('\rho_{gas}/\rho_{vir}','Fontsize',12);
%     title(sprintf('Gas Density, %s{vir}=%s M_{sun}/Mpc^3', texlabel('rho'),num2str(rhovir,'%1.2d')),'Fontsize',12);
%     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);
%     
%     %subplot(2,2,4);
%     %semilogx(r_prof/rvir,vrprof./vvir);grid;xlim([1e-3 3]);%ylim([ ]);
%     %xlabel('r/R_{vir}','Fontsize',12);
%     %ylabel('V_r/V_{vir}','Fontsize',12);
%     %title(sprintf('Radial Velocity, V_{vir}=%s km/sec', num2str(vvir,'%1.2d')),'Fontsize',12);
%     %set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);  
%     
%     subplot(2,2,4);
%     loglog(r_prof/rvir,s_prof./(tvir./rhovir.^(2/3)) );grid;xlim([1e-3 3]);ylim([1e-4 1e2]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('S/S_{vir}','Fontsize',12);
%     title(sprintf('Entropy, S_{vir}=%s K Mpc^2 /{M_{sun}}^{2/3}',num2str((tvir./rhovir.^(2/3)),'%1.2d')),'Fontsize',12);
%     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);
%     
%     if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_all_profs.png',result_dir,clustername));
%     end
%     
%     end
%         
%     close all
% %end
% %save('mat_files/stk_profs.mat','stack');
% 
%   
