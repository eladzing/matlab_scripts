%% load 

cluster='CL3';
xtxt=-1.5;
ytxt1=0.07;ytxt2=0.06;


figure;

load(sprintf('%s/%s_tfluxbird_phd_clean_vzone2.mat','mat_files',cluster));
%load(sprintf('%s/%s_tfluxbird_phd_outer_clean_vzone2.mat','mat_files',cluster));

%% plot bird
subplot(1,2,1);

cnorm=squeeze(tro(:,:,1));
imagesc(tlim,rolim,log10(tro(:,:,1)./sum(cnorm(:))));
%bar=colorbar;
%set(get(bar,'Title'),'String',sprintf('$\mathrm{log}\left(\frac{M}{M_{zone}}\right)\,;\,pixel=%s$',num2str(binorm,'%1.2g')),'Fontsize',12,'Interpreter','latex')
%set(get(bar,'Title'),'String','$\mathrm{log}\left(M/M_{zone}\right)$','Fontsize',12,'Interpreter','latex')
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
%caxis([-6.2 -2.0])
%set(gca,'Xdir','normal');
set(gca,'Ydir','normal');
xlabel('$\mathrm{log}\left(T/T_{\mathrm{vir}}\right)$','Fontsize',12,'Interpreter','latex');
ylabel('$\dot M/M_{gas}(R_{\mathrm{vir}})/d\Omega\,[1/\mathrm{Gyr}]$','Fontsize',12,'Interpreter','latex');grid;

text(xtxt,ytxt1,'Inner Zone','fontsize',12,'interpreter','latex')
text(xtxt,ytxt2,'$0.01<r/R_{\mathrm{vir}}<0.2$','fontsize',12,'interpreter','latex')


%title(sprintf( '%s Mass Histogram (%s<r/R_{vir}<%s)',CLUSTER,num2str(rmfac,2),num2str(rxfac,2)),'Fontsize',12);
set(gca,...
    'Position',[0.1 0.15 0.37 0.7],...
    'FontSize',12);
    %'PlotBoxAspectRatio',[1 1 1],...
    %DataAspectRatio',[1 1 1]),...
    %'FontSize',12);
    
  
subplot(1,2,2);

load(sprintf('%s/%s_tfluxbird_phd_outer_clean_vzone2.mat','mat_files',cluster));
cnorm=squeeze(tro(:,:,1));
imagesc(tlim,rolim,log10(tro(:,:,1)./sum(cnorm(:))));

%set(get(bar,'Title'),'String',sprintf('$\mathrm{log}\left(\frac{M}{M_{zone}}\right)\,;\,pixel=%s$',num2str(binorm,'%1.2g')),'Fontsize',12,'Interpreter','latex')
%set(get(bar,'Title'),'String','$\mathrm{log}\left(M/M_{zone}\right)$','Fontsize',12,'Interpreter','latex')
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
caxis([-6.2 -2.0])
%set(gca,'Xdir','normal');
set(gca,'Ydir','normal');
grid;
xlabel('$\mathrm{log}\left(T/T_{\mathrm{vir}}\right)$','Fontsize',12,'Interpreter','latex');
%ylabel('$\dot M/M_{gas}(R_{\mathrm{vir}})/d\Omega\,[1/\mathrm{Gyr}]$','Fontsize',12,'Interpreter','latex');grid;
%title(sprintf( '%s Mass Histogram (%s<r/R_{vir}<%s)',CLUSTER,num2str(rmfac,2),num2str(rxfac,2)),'Fontsize',12);
set(gca,...
    'Position',[0.53 0.15 0.37 0.7],...
    'FontSize',12);
    %'PlotBoxAspectRatio',[1 1 1],...
    %'DataAspectRatio',[1 1 1]);
 bar=colorbar('position',[0.91 0.15 0.04 0.7]);
 %set(get(bar,'Title'),'String',sprintf('$\mathrm{log}\left(\frac{M}{M_{zone}}\right)\,;\,pixel=%s$',num2str(binorm,'%1.2g')),'Fontsize',12,'Interpreter','latex')
set(get(bar,'Title'),'String','$\mathrm{log}\left(M/M_{zone}\right)\quad$','Fontsize',12,'Interpreter','latex');
text(xtxt,ytxt1,'Outer Zone','fontsize',12,'interpreter','latex')
text(xtxt,ytxt2,'$0.2<r/R_{\mathrm{vir}}<1.2$','fontsize',12,'interpreter','latex')




