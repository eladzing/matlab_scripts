%list1=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%list1=[101 102 103 104 105 106 107 3 6 7 9 10 11 14 24];

%list2=[103 105 106 107];

result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout';


load('mat_files/stk_penen.mat','stackpen');

plotflag='noplot';
pflag='noprint';
 
%vzone=[0.02 0.2];
r_dens=0.01;
r_md=0.2;

ropar=[];
flxpar=[];
cname=[];

%     stackpen{id,1}=clustername;
%     stackpen{id,2}=rp;
%     stackpen{id,3}=flux1;
%     stackpen{id,4}=flux2;
%     stackpen{id,5}=flux3;
%     stackpen{id,6}=flux4;
%     stackpen{id,7}=rps;
%     stackpen{id,8}=[mvir,rvir,Mg200,rhovir,vvir,fluxnorm];

for id=1:16 %length(stackpen{:,1} )
    
    cname{end+1}=stackpen{id,1};
    mvir=stackpen{id,9}(1);
    rvir=stackpen{id,9}(2);
    Mg200=stackpen{id,9}(3);
    rhovir=stackpen{id,9}(4);
    vvir=stackpen{id,9}(5);
    %flnorm=stackpen{id,9}(6);
        
    
    
    rop=stackpen{id,8};
    rp=stackpen{id,2};
    fluxnorm=0.056.*(mvir./1e13).^0.15;
    flx1=stackpen{id,3};
    flx2=stackpen{id,4};
    flx3=stackpen{id,5};
    
    flxpar(end+1,1)=interp1(rp,flx1,r_md,'linear');
    flxpar(end+1,2)=interp1(rp,flx2,r_md,'linear');    
    flxpar(end+1,3)=interp1(rp,flx3,r_md,'linear');
    
    ropar(end+1)=interp1(rp,rop,r_dens,'linear');
    
end
    
    
    
    
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
    % text(0.1,0.7,sprintf('M_{gas}=%s',num2str(Mg200,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
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
%      text(0.1,0.7,sprintf('M_{gas}=%s M_{sun}',num2str(Mg200,'%1.2d')) , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Fontsize',12);
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
