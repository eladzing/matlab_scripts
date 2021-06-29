
load mat_files/stk_profs_a1.mat

flist=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];

%      1   2   3   4   5   6   7  8 9 10111213 14 15 16   
pflag='noprint';
cflist=[104 105 107 3 6 7 10 14;4 5 7 8 10 11 13 15];
ncflist=[101 102 103 106 5 9 11 24;1 2 3 6 9 12 14 16];
result_dir='~/work/cold_flows/datacube/printout/';

%  stack{id,2}=r_prof/rvir; %r
%     stack{id,3}=m_dot./Mg200; %Flux profile
%     stack{id,4}=t_prof./tvir; %temperature profile 
%     stack{id,5}=roprof./rhovir; %density profile 
%     stack{id,6}=vrprof./vvir; %radial velocity profile 
%     stack{id,7}=s_prof./(tvir./rhovir.^(2/3)); %full mass flux
%     stack{id,8}=[mvir,rvir,Mg200,tvir,rhovir,vvir,(tvir./rhovir.^(2/3))];

rp=0.001:0.001:2.6;

for jj=1:1
    switch jj
        case 1
        lst=1:1:16; 
        printag='all'
        titletag='All';
%         case 2;
%         lst=cflist(2,:);
%         printag='cf'
%         titletag='CF';
%         case 3
%         lst=ncflist(2,:);
%         printag='ncf'
%         titletag='NCF';
    end
    n=0;  
    mdf=[];
    tp=[];
    rop=[];
    vrp=[];
    sp=[];
    for i=lst %%1:length(lst);
        %clname=sprintf('CL%d',cflist(i));
        %for j=1:size(stack,1)
        %ind=find(strcmp(stack{:,1}
        %if i~=9
        normv=stack{i,8};
        n=n+1;
        clname{n}=sprintf('CL%s',num2str(flist(i)));
        mdf(end+1,:)=interp1(stack{i,2},stack{i,3},rp,'linear').*1e9;   %%./(0.056.*(normv(1)./1e13).^0.15).*1e9;
        tp(end+1,:)=interp1(stack{i,2},stack{i,4},rp,'linear');
        rop(end+1,:)=interp1(stack{i,2},stack{i,5},rp,'linear');
        vrp(end+1,:)=interp1(stack{i,2},stack{i,6},rp,'linear');
        sp(end+1,:)=interp1(stack{i,2},stack{i,7},rp,'linear');
        %end
    end

%%plot stacks
    figure;
    
    
    subplot(4,2,1);
    semilogx(rp,mdf','linewidth',1);grid;xlim([1e-2 3]);ylim([-0.38 0.2]);
    line([0.2 0.2],[-1e5 1e5],'Color','g','LineWidth',1.2); 
    %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
    ylabel('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$','Fontsize',12,'Interpreter','latex');% / M_{vir} [Gyr^{-1}] $}','Fontsize',12,'Interpreter','latex');
    %title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
    line([1e-3 1e1],[0 0],'Color','k');
    
    %legend(clname,'location','bestoutside');
    set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'xticklabel','','Fontsize',12,'linewidth',1,...
        'Position',[0.15 0.77 0.42 0.21]);
        
    fl1=mean(mdf,1);%fl2=sum(mdf2,1);
    v1=std(mdf,0,1);%v2=std(mdf2,0,1);
   
    subplot(4,2,2);
    semilogx(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([-0.38 0.2]);hold on;
    semilogx(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
    hold off;

    % xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
    % ylabel('$\dot{M}/M_{\mathrm{vir}} [\mathrm{Gyr}^{-1}]$','Fontsize',12);
     
     line([1e-3 1e1],[0 0],'Color','k');
     line([0.2 0.2],[-1e5 1e5],'Color','g','LineWidth',1.2);   
     set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.4 -0.2 0 0.2],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
         'Position',[0.57 0.77 0.42 0.21]);
     %laprint(gcf,'flx_stk','options','factory','width',9,'factor',1,'scalefonts','on');
    
    
    %figure;
    subplot(4,2,3);
    loglog(rp,tp','linewidth',1);grid;xlim([1e-2 3]);ylim([5e-2 10]);
    line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
    %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
    ylabel('$T/T_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
    set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'xticklabel','','Fontsize',12,'linewidth',1,...
        'Position',[0.15 0.56 0.42 0.21]);
    %laprint(gcf,'tee_all','options','factory','width',9,'factor',1,'scalefonts','on');
    
    
    fl1=mean(tp,1);%fl2=sum(mdf2,1);
    v1=std(tp,0,1);%v2=std(mdf2,0,1);
   
    %figure;
    subplot(4,2,4);
    loglog(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([5e-2 10]);hold on;
    loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
    line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
    hold off;
    %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
    %ylabel('$T/T_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
    set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
        'Position',[0.57 0.56 0.42 0.21]);
    %laprint(gcf,'tee_stk','options','factory','width',9,'factor',1,'scalefonts','on');

    %figure;
    subplot(4,2,5)
    loglog(rp,rop','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-3 1e3]);
    line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
    %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
    ylabel('$\rho_{gas}/\rho_{\mathrm{mean}}$','Fontsize',12,'Interpreter','latex');   
    set(gca,'XTick', [1e-2 1e-1 1 5],'YTick',[1e-2 1 1e2 ],'xticklabel','','Fontsize',12,'linewidth',1,...
       'Position',[0.15 0.35 0.42 0.21] );
    %laprint(gcf,'rho_all','options','factory','width',9,'factor',1,'scalefonts','on');
    
    fl1=mean(rop,1);%fl2=sum(mdf2,1);
    v1=std(rop,0,1);%v2=std(mdf2,0,1);
   
    %figure;
    subplot(4,2,6);
    loglog(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-3 1e3]);hold on;
    loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
    line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
    hold off;
    %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12);
    %ylabel('$\rho_{gas}/\rho_{\mathrm{vir}}$','Fontsize',12);
    set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-2 1 1e2],'xticklabel','','yticklabel','','Fontsize',12,'linewidth',1,...
     'Position',[0.57 0.35 0.42 0.21]);
    %laprint(gcf,'rho_stk','options','factory','width',9,'factor',1,'scalefonts','on');

   %figure;
    subplot(4,2,7);
    loglog(rp,sp','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-2 1e1]);
   line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
    xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');  
   ylabel('$S/S_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');   
   set(gca,'XTick', [1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'Fontsize',12,'linewidth',1,...
   'Position',[0.15 0.14 0.42 0.21]);
   % laprint(gcf,'ent_all','options','factory','width',9,'factor',1,'scalefonts','on');
    
    
    fl1=mean(sp,1);%fl2=sum(mdf2,1);
    v1=std(sp,0,1);%v2=std(mdf2,0,1);
   
    %figure;
    subplot(4,2,8);
    loglog(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-2 1e1]);hold on;
    loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
    line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
    hold off;
      xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
    %ylabel('$S/S_{\mathrm{vir}}$','Fontsize',12);
   
     set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'yticklabel','','Fontsize',12,'linewidth',1,...
         'Position',[0.57 0.14 0.42 0.21]);
     %laprint(gcf,'ent_stk','options','factory','width',9,'factor',1,'scalefonts','on');
    
     if strcmp(pflag,'print')
         %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
         saveas(gcf,sprintf('%s/%s_stack_profs.png',result_dir,'multi'));
     end
     
    clear normv clname mdf tp rop vrp sp;
end
    
   
    





    
