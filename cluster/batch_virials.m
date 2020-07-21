list=[101 102 103 104 105 106 107  3  5  6  7  9 10 11 14 24];
%mask=[ 0   0   0   0   0   0   0   1  1  1  1  1  1  1  1  1];
mask=ones(size(list));
cluster_virials=cell(length(list),4);

aexp='a06';

outputdir='printout';%../cluster_printout';

virs=[];
virs_200=[];
virs_500=[];
virs_old=[];

explain='contains 1- cluster name 2-4:[rv mv vv tv deltavir] for vir,200 & 500';

for i=1:length(list);
    if~mask(i) 
        continue;
    end
    new_env(sprintf('CL%d',list(i)),'csf',aexp);
    
   global CLUSTER
   
   [r m v t d]=calc_virials();
   
   virs(end+1,1:5)=[r m v t d];
  [r m v t d] =calc_virials(200);
   virs_200(end+1,1:5)=[r m v t d];
  [r m v t d]=calc_virials(500);
  virs_500(end+1,1:5)=[r m v t d];
 
   
   %virs_old(end+1,:)=[get_rvir() get_mvir() get_vvir() get_tvir()];
  
   cluster_virials{i,1}=CLUSTER;
   cluster_virials{i,2}=virs(end,:);
   cluster_virials{i,3}=virs_200(end,:);
   cluster_virials{i,4}=virs_500(end,:);
   
   
end
cluster_virials{end+1,1}=explain;
ii=1:sum(mask);
labcell={'101','102','103','104','105','106','107','3','5','6','7','9','10','11','14','24'};

% figure
% 
% 
% subplot(3,1,1:2)
% plot(ii,virs_old(:,1),'r*',ii,virs(:,1),'bo','markersize',8)
% ylabelmine('$R_{vir}\,[\mathrm{Mpc}]$');
% format_ticks(gca,labcell,[],ii);
% titlemine('New vs. Old virial Calculation')
% 
% subplot(3,1,3)
% plot(ii,virs(:,1)./virs_old(:,1)-1,'bo','markersize',8)
% ylabelmine('residuals');
% format_ticks(gca,labcell,[],ii);
% xlabelmine('Halo');

%saveas(gcf,sprintf('%s/rvir_comparison_%s.png',outputdir,aexp));

% figure
% 
% 
% subplot(3,1,1:2)
% plot(ii,virs_old(:,2),'r*',ii,virs(:,2),'bo','markersize',8)
% ylabelmine('$M_{vir}\, [\mathrm{M_\odot}]$');
% format_ticks(gca,labcell,[],ii);
% titlemine('New vs. Old virial Calculation')
% 
% subplot(3,1,3)
% plot(ii,virs(:,1)./virs_old(:,1)-1,'bo','markersize',8)
% 
% format_ticks(gca,labcell,[],ii);
% 
% xlabelmine('Halo');
% ylabelmine('residuals');
% 
% %saveas(gcf,sprintf('%s/mvir_comparison_%s.png',outputdir,aexp));
% 
% 
% figure
% 
% subplot(3,1,1:2)
% plot(ii,virs_old(:,3),'r*',ii,virs(:,3),'bo','markersize',8)
% format_ticks(gca,labcell,[],ii);
% ylabelmine('$V_{vir}\,[\mathrm{km/sec}]$');
% titlemine('New vs. Old virial Calculation')
% 
% subplot(3,1,3)
% plot(ii,virs(:,1)./virs_old(:,1)-1,'bo','markersize',8)
% 
% format_ticks(gca,labcell,[],ii);
% 
% xlabelmine('Halo');
% ylabelmine('residuals');
% 
% %saveas(gcf,sprintf('%s/vvir_comparison_%s.png',outputdir,aexp));
% 
% figure
% 
% subplot(3,1,1:2)
% plot(ii,virs_old(:,4),'r*',ii,virs(:,4),'bo','markersize',8)
% ylabelmine('$T_{vir}\,[\mathrm{K}]$');
% format_ticks(gca,labcell,[],ii);
% titlemine('New vs. Old virial Calculation')
% 
% subplot(3,1,3)
% plot(ii,virs(:,1)./virs_old(:,1)-1,'bo','markersize',8)
% 
% format_ticks(gca,labcell,[],ii);
% xlabelmine('Halo');
% ylabelmine('residuals');
% 
%saveas(gcf,sprintf('%s/tvir_comparison_%s.png',outputdir,aexp));

figure
plot(ii,virs(:,1),'bo',ii,virs_200(:,1),'ro',ii,virs_500(:,1),'go','markersize',7)
format_ticks(gca,labcell,[],ii);
xlabelmine('Halo');
ylabelmine('$R_{vir}\,[\mathrm{Mpc}]$');
h=legend('$\Delta_{vir}=337$','$\Delta_{vir}=200$','$\Delta_{vir}=500$');
set(h,'Interpreter','latex','Fontsize',12);
titlemine('$R_{vir}$ for varying $\Delta_{vir}$')

%saveas(gcf,sprintf('%s/rvir_%s.png',outputdir,aexp));

figure
plot(ii,virs(:,2),'bo',ii,virs_200(:,2),'ro',ii,virs_500(:,2),'go','markersize',7)
format_ticks(gca,labcell,[],ii);
xlabelmine('Halo');
ylabelmine('$M_{vir}\,[\mathrm{M\odot}]$');
h=legend('$\Delta_{vir}=337$','$\Delta_{vir}=200$','$\Delta_{vir}=500$');
set(h,'Interpreter','latex','Fontsize',12);
titlemine('$M_{vir}$ for varying $\Delta_{vir}$')

%saveas(gcf,sprintf('%s/mvir_%s.png',outputdir,aexp));

 figure
plot(ii,virs(:,3),'bo',ii,virs_200(:,3),'ro',ii,virs_500(:,3),'go','markersize',7)
format_ticks(gca,labcell,[],ii);
xlabelmine('Halo');
ylabelmine('$V_{vir}\,[\mathrm{km/sec}]$');
h=legend('$\Delta_{vir}=337$','$\Delta_{vir}=200$','$\Delta_{vir}=500$');
set(h,'Interpreter','latex','Fontsize',12);
titlemine('$V_{vir}$ for varying $\Delta_{vir}$')

%saveas(gcf,sprintf('%s/vvir_%s.png',outputdir,aexp));

figure
plot(ii,virs(:,4),'bo',ii,virs_200(:,4),'ro',ii,virs_500(:,4),'go','markersize',7)
format_ticks(gca,labcell,[],ii);
xlabelmine('Halo');
ylabelmine('$T_{vir}\,[\mathrm{K}]$');
h=legend('$\Delta_{vir}=337$','$\Delta_{vir}=200$','$\Delta_{vir}=500$');
set(h,'Interpreter','latex','Fontsize',12);
titlemine('$T_{vir}$ for varying $\Delta_{vir}$')
   
%saveas(gcf,sprintf('%s/tvir_%s.png',outputdir,aexp));  
   
   
