
load mat_files/stk_profs_a1.mat

flist=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
%      1    2   3   4   5   6   7  8 9 10111213 14 15 16   
pflag='print';

% cflist=[104 105 107 3 6 7 10 14;4 5 7 8 10 11 13 15];
% ncflist=[101 102 103 106 5 9 11 24;1 2 3 6 9 12 14 16];

% r_flist=[104 3 5 7 10 14;4 8 9 11 13 15]; %relaxed
% ur_flist=[101 102 103 105 106 107 6 9 11 24;1 2 3 5 6 7 10 12 14 16]; %unrelaxed
% 
% d_flist=[101 103 105 107 3 6 14;1 3 5 7 8 10 15]; %deep     
% s_flist=[102 104 5 7 10 24;2 4 9 11 13 16]; %shallow
% 
% b_flist=[102 103 106 5 9 24;2 3 6 9 12 16]; %bump
% nb_flist=[101 104 105 107 3 6 7 10 11 14;1 4 5 7 8 10 11 13 14 15]; %no bump

result_dir='~/work/cold_flows/datacube/printout/';

%  stack{id,2}=r_prof/rvir; %r
%     stack{id,3}=m_dot./Mg200; %Flux profile
%     stack{id,4}=t_prof./tvir; %temperature profile 
%     stack{id,5}=roprof./rhovir; %density profile 
%     stack{id,6}=vrprof./vvir; %radial velocity profile 
%     stack{id,7}=s_prof./(tvir./rhovir.^(2/3)); %full mass flux
%     stack{id,8}=[mvir,rvir,Mg200,tvir,rhovir,vvir,(tvir./rhovir.^(2/3))];

rp=0.001:0.001:2.6;


    n=0;  
    mdf=[];
    tp=[];
    rop=[];
    vrp=[];
    sp=[];
    for i=lst %%1:length(lst);
        %clname=sprintf('CL%d',flist(i));
        %for j=1:size(stack,1)
        %ind=find(strcmp(stack{:,1}
        %if i~=14
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
    figure('PaperType','A4','PaperOrientation','landscape','PaperPosition',[0.5 0.5 10.5 8]);
    %title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
    for j=lst
        len=0.22;
        a=mod(j-1,4);
        b=4-ceil(j/4);
        left=0.1+len.*a ;
        bottom= 0.1+len.*b;
        subplot('Position',[left bottom (len-0.02) (len-0.02)])
        %subplot(4,4,j);
        semilogx(rp,mdf(j,:),'linewidth',1);grid;xlim([1e-2 3]);ylim([-0.38 0.2]);
        %set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1);
        line([1e-3 1e1],[0 0],'Color','k');
        text(0.021,0.11,sprintf('%s',clname{j}),'Fontsize',9);
        switch (a)
            case 0
                %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',10,'Interpreter','latex');
                ylabel('$\dot{M}/M_{g} [\mathrm{Gyr}^{-1}]$','Fontsize',10,'Interpreter','latex');
                set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1,...
                    'xticklabel','');
            otherwise
                set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1,...
                    'yticklabel','','xticklabel','');
        end
        switch (b)
            case 0
                xlabel('$log(r/R_{\mathrm{vir}})$','Fontsize',10,'Interpreter','latex');
                %ylabel('$\dot{M}/M_{g} [\mathrm{Gyr}^{-1}]$','Fontsize',10,'Interpreter','latex');
                set(gca,'xticklabel',[0 -2 -1]);
            otherwise
                %set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1,...
                %    'yticklabel','','xticklabel','');
        end
    end
    
    if strcmp(pflag,'print')
%          %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
          saveas(gcf,sprintf('%s/%s_tile_profs.png',result_dir,'mflx'));
    end
    
%     
%     
%     
%     subplot(2,1,1);
%     semilogx(rp,mdf','linewidth',1);grid;xlim([1e-2 3]);ylim([-0.38 0.2]);
%     line([0.2 0.2],[-1e5 1e5],'Color','g','LineWidth',1.2); 
%     %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex','Interpreter','latex');
%     ylabel('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$','Fontsize',12,'Interpreter','latex');% / M_{vir} [Gyr^{-1}] $}','Fontsize',12,'Interpreter','latex');
%     title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
%     line([1e-3 1e1],[0 0],'Color','k');
%     
%     legend(clname,'location','bestoutside');
%     
%     set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'xticklabel','','Fontsize',12,'linewidth',1,...
%         'Position',[0.12 0.52 0.7 0.40]);
%         
%     fl1=mean(mdf,1);%fl2=sum(mdf2,1);
%     v1=std(mdf,0,1);%v2=std(mdf2,0,1);
%     %figure;
%     subplot(2,1,2);
%     semilogx(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([-0.38 0.2]);hold on;
%     semilogx(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
%     hold off;
% 
%      xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%      ylabel('$\dot{M}/M_{gas} [1/\mathrm{Gyr}]$','Fontsize',12,'Interpreter','latex');
%      %%title(sprintf('%s Clusters Mass Flux',titletag),'Fontsize',12);
%      line([1e-3 1e1],[0 0],'Color','k');
%      line([0.2 0.2],[-1e5 1e5],'Color','g','LineWidth',1.2);   
%      set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-0.3 -0.2 -0.1 0 0.1 0.2],'Fontsize',12,'linewidth',1,...
%          'Position',[0.12 0.10 0.7 0.40]);
%      if strcmp(pflag,'print')
%          %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
%          saveas(gcf,sprintf('%s/%s_stack_profs_%s.png',result_dir,'mflx',printag));
%      end
%     
%     figure;
%     subplot(2,1,1);
%     semilogx(rp,vrp','linewidth',1);grid;xlim([1e-2 3]);ylim([-1 0.5]);
%     line([0.2 0.2],[-1e5 1e5],'Color','g','LineWidth',1.2); 
%     %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     ylabel('$v_r/v_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');% / M_{vir} [Gyr^{-1}] $}','Fontsize',12,'Interpreter','latex');
%     title(sprintf('%s Clusters Radial Velocity',titletag),'Fontsize',12);
%     line([1e-3 1e1],[0 0],'Color','k');
%     legend(clname,'location','bestoutside');
%     
%     set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-1 -0.5 0 0.5],'Fontsize',12,'linewidth',1,...
%         'Position',[0.12 0.52 0.7 0.40]);
%         
%     fl1=mean(vrp,1);%fl2=sum(mdf2,1);
%     v1=std(vrp,0,1);%v2=std(mdf2,0,1);
%    
%     subplot(2,1,2);
%     semilogx(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([-1 0.5]);hold on;
%     semilogx(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
%     hold off;
% 
%      xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%      ylabel('$v_r/v_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%      %title(sprintf('%s Clusters Radial Velocity',titletag),'Fontsize',12);
%      line([1e-3 1e1],[0 0],'Color','k');
%      line([0.2 0.2],[-1e5 1e5],'Color','g','LineWidth',1.2);   
%      set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[-1 -0.5 0 0.5],'Fontsize',12,'linewidth',1,...
%          'Position',[0.12 0.10 0.7 0.40]);
%      
%      if strcmp(pflag,'print')
%          %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
%          saveas(gcf,sprintf('%s/%s_stack_profs_%s.png',result_dir,'vr',printag));
%      end
%      
%     
%     figure;
%     subplot(2,1,1);
%     loglog(rp,tp','linewidth',1);grid;xlim([1e-2 3]);ylim([5e-2 10]);
%     line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
%     %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     ylabel('$T/T_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     title(sprintf('%s Clusters Temperature',titletag),'Fontsize',12);
%     legend(clname,'location','bestoutside');
%     set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'xticklabel','','Fontsize',12,'linewidth',1,...
%         'Position',[0.12 0.52 0.7 0.40]);
%     %laprint(gcf,'tee_all','options','factory','width',9,'factor',1,'scalefonts','on');
%     
%     
%     fl1=mean(tp,1);%fl2=sum(mdf2,1);
%     v1=std(tp,0,1);%v2=std(mdf2,0,1);
%    
%      subplot(2,1,2);
%     loglog(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([5e-2 10]);hold on;
%     loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
%     line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
%     hold off;
%     xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     ylabel('$T/T_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     %title(sprintf('%s Clusters Temperature',titletag),'Fontsize',12);
%     set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'Fontsize',12,'linewidth',1,...
%         'Position',[0.12 0.10 0.7 0.40]);
%     if strcmp(pflag,'print')
%          %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
%          saveas(gcf,sprintf('%s/%s_stack_profs_%s.png',result_dir,'tmp',printag));
%      end
%      
% 
%     figure;
%     subplot(2,1,1);
%     loglog(rp,rop','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-3 1e3]);
%     line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
%     %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     ylabel('$\rho_{gas}/\rho_{\mathrm{mean}}$','Fontsize',12,'Interpreter','latex');   
%     title(sprintf('%s Clusters Density',titletag),'Fontsize',12);
%     legend(clname,'location','bestoutside');
%     set(gca,'XTick', [1e-2 1e-1 1 5],'YTick',[1e-2 1 1e2 ],'xticklabel','','Fontsize',12,'linewidth',1,...
%        'Position',[0.12 0.52 0.7 0.40]);
%     
%     
%     
%     fl1=mean(rop,1);%fl2=sum(mdf2,1);
%     v1=std(rop,0,1);%v2=std(mdf2,0,1);
%    
%     subplot(2,1,2)
%     loglog(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-3 1e3]);hold on;
%     loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
%     line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
%     hold off;
%     xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     ylabel('$\rho_{gas}/\rho_{\mathrm{mean}}$','Fontsize',12,'Interpreter','latex');  
%     %title(sprintf('%s Clusters Density',titletag),'Fontsize',12);
%     set(gca,'XTick',[1e-3 1e-2 1e-1 1 5],'YTick',[1e-2 1 1e2],'Fontsize',12,'linewidth',1,...
%      'Position',[0.12 0.10 0.7 0.40]);
%     
%     if strcmp(pflag,'print')
%          %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
%          saveas(gcf,sprintf('%s/%s_stack_profs_%s.png',result_dir,'ro',printag));
%      end
% 
%    figure;
%    subplot(2,1,1);
% 
%     %subplot(4,2,7);
%     loglog(rp,sp','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-2 1e1]);
%    line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
%     %xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');  
%    ylabel('$S/S_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');   
%    title(sprintf('%s Clusters Entropy',titletag),'Fontsize',12);
%    legend(clname,'location','bestoutside');
%    set(gca,'XTick', [1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'xticklabel','','Fontsize',12,'linewidth',1,...
%    'Position',[0.12 0.52 0.7 0.40]);
%    
%     fl1=mean(sp,1);%fl2=sum(mdf2,1);
%     v1=std(sp,0,1);%v2=std(mdf2,0,1);
%    
%     subplot(2,1,2);
%     loglog(rp,fl1,'-b','linewidth',1);grid;xlim([1e-2 3]);ylim([1e-2 1e1]);hold on;
%     loglog(rp,[fl1-v1./2;fl1+v1./2],'-.b','linewidth',1);%grid;xlim([1e-3 10]);ylim([-3 1]);
%     line([0.2 0.2],[1e-5 1e5],'Color','g','LineWidth',1.2);  
%     hold off;
%       xlabel('$r/R_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%     ylabel('$S/S_{\mathrm{vir}}$','Fontsize',12,'Interpreter','latex');
%    %title(sprintf('%s Clusters Entropy',titletag),'Fontsize',12);
%      set(gca,'XTick', [1e-3 1e-2 1e-1 1 5],'YTick',[1e-1 1 1e1],'Fontsize',12,'linewidth',1,...
%          'Position',[0.12 0.10 0.7 0.40]);
%      
%     if strcmp(pflag,'print')
%          %saveas(gcf,sprintf('%s/%s_stack_profs',result_dir,'multi'),'eps');
%          saveas(gcf,sprintf('%s/%s_stack_profs_%s.png',result_dir,'ent',printag));
%      end
%           
%     clear normv clname mdf tp rop vrp sp;
%     %close all

    
   
    





    
