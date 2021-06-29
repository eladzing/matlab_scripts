
load mat_files/stk_profs.mat

flist=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%      1   2   3   4   5   6   7  8 9 10111213 14 15 16   
pflag='noprint';
cflist=[104 105 107 6 10 14;4 5 7 10 13 15];
ncflist=[101 102 103 106 3 5 7 9 11 24;1 2 3 6 8 9 11 12 14 16];
result_dir='/home/carina/eladzing_a/cold_flows/tit3/printout';

%  stack{id,2}=r_prof/rvir; %r
%     stack{id,3}=m_dot./Mg200; %Flux profile
%     stack{id,4}=t_prof./tvir; %temperature profile 
%     stack{id,5}=roprof./rhovir; %density profile 
%     stack{id,6}=vrprof./vvir; %radial velocity profile 
%     stack{id,7}=s_prof./(tvir./rhovir.^(2/3)); %full mass flux
%     stack{id,8}=[mvir,rvir,Mg200,tvir,rhovir,vvir,(tvir./rhovir.^(2/3))];

rp=0.001:0.001:2.5;

for jj=1:3
    switch jj
        case 1
        lst=1:1:16; 
        printag='all'
        titletag='All';
        case 2;
        lst=cflist(2,:);
        printag='cf'
        titletag='CF';
        case 3
        lst=ncflist(2,:);
        printag='ncf'
        titletag='NCF';
    end
    n=0;  
    mdf1=[];
    mdf2=[];
    tp=[];
    rop=[];
    vrp=[];
    sp=[];
    for i=lst %%1:length(lst);
        %clname=sprintf('CL%d',cflist(i));
        %for j=1:size(stack,1)
        %ind=find(strcmp(stack{:,1}
        normv=stack{i,8};
        n=n+1;
        clname{n}=sprintf('CL%s',num2str(flist(i)));
        mdf1(end+1,:)=interp1(stack{i,2},stack{i,3},rp,'linear')./(0.056.*(normv(1)./1e13).^0.15).*1e9;
        %mdf2(end+1,:)=interp1(stack{i,2},stack{i,3},rp,'linear');
        tp(end+1,:)=interp1(stack{i,2},stack{i,4},rp,'linear');
        rop(end+1,:)=interp1(stack{i,2},stack{i,5},rp,'linear');
        vrp(end+1,:)=interp1(stack{i,2},stack{i,6},rp,'linear');
        sp(end+1,:)=interp1(stack{i,2},stack{i,7},rp,'linear');
    end

%% add stacks


fl1=mean(mdf1,1);%fl2=sum(mdf2,1);
v1=std(mdf1,0,1);v2=std(mdf2,0,1);

figure;

semilogx(rp,fl1,'-b');grid;xlim([1e-3 10]);ylim([-3 1]);hold on;
semilogx(rp,[fl1-v1./2;fl1+v1./2],'-.b');%grid;xlim([1e-3 10]);ylim([-3 1]);
hold off;
%semilogx(rp,fl1,'-b');grid;xlim([1e-3 10]);%ylim([-6 2]);
     xlabel('r/R_{vir}','Fontsize',12);
     ylabel('dM/dt','Fontsize',12);
     title(sprintf('%s Clusters Stacked Flux',titletag),'Fontsize',12);
     line([1e-3 1e1],[0 0],'Color','k');
%     legend(clname,'location','bestoutside');
     set(gca,'Fontsize',12);
    
   %clear normv clname mdf tp rop vrp sp fl1;

if strcmp(pflag,'print')
         saveas(gcf,sprintf('%s/%s_addstack_flux_prof.png',result_dir,printag));
      end
end

% 
% 
% %%plot stacks
%     figure;
%     
%     semilogx(rp,mdf1');grid;xlim([1e-3 10]);ylim([-6 2]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('dM/dt [virial dm/dt]','Fontsize',12);
%     title(sprintf('%s Clusters Stacked Flux',titletag),'Fontsize',12);
%     line([1e-3 1e1],[0 0],'Color','k');
%     legend(clname,'location','bestoutside');
%     set(gca,'Fontsize',12);
%    
%     
%      if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_stack_flux_prof.png',result_dir,printag));
%      end
%         
%      figure;
%     
%     semilogx(rp,md1f');grid;xlim([1e-3 10]);ylim([-6 2]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('dM/dt [virial dm/dt]','Fontsize',12);
%     title(sprintf('%s Clusters Stacked Flux',titletag),'Fontsize',12);
%     line([1e-3 1e1],[0 0],'Color','k');
%     legend(clname,'location','bestoutside');
%     set(gca,'Fontsize',12);
%    
%     
%      if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_stack_flux_prof.png',result_dir,printag));
%      end
%     
%      
%      
%      
%     figure;
%     semilogx(rp,vrp');grid;xlim([1e-3 3]);ylim([-1 0.25]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel(sprintf('V_r/V_{vir}'),'Fontsize',12);
%     title(sprintf('%s Clusters Stacked Radial Velocity',titletag),'Fontsize',12)
%     legend(clname,'location','bestoutside');
%     set(gca,'Fontsize',12);  
%     if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_stack_vr_prof.png',result_dir,printag));
%     end
%     
%         
%     
%     figure;
%     loglog(rp,tp');grid;xlim([1e-3 3]);ylim([5e-2 6]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('T/T_{vir}','Fontsize',12);
%     title(sprintf('%s Clusters Stacked Temperature',titletag),'Fontsize',12);
%     legend(clname,'location','bestoutside');
%     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);
%     if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_stack_t_prof.png',result_dir,printag));
%     end
%     
%     figure;
%     loglog(rp,rop');grid;xlim([1e-3 3]);ylim([8e-3 1e4]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('\rho_{gas}/\rho_{vir}','Fontsize',12);
%     title(sprintf('%s Clusters Stacked Density',titletag),'Fontsize',12);
%     legend(clname,'location','bestoutside');
%     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);
%     if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_stack_ro_prof.png',result_dir,printag));
%     end
%     
%    figure;
%      loglog(rp,sp');grid;xlim([1e-3 3]);ylim([1e-4 1e2]);
%     xlabel('r/R_{vir}','Fontsize',12);
%     ylabel('S/S_{vir}','Fontsize',12);
%     title(sprintf('%s Clusters Stacked  Entropy',titletag),'Fontsize',12);
%     legend(clname,'location','bestoutside');
%     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'Fontsize',12);
%     if strcmp(pflag,'print')
%         saveas(gcf,sprintf('%s/%s_stack_s_prof.png',result_dir,printag));
%     end
%     clear normv clname mdf tp rop vrp sp;
% end
    
   
    





    
