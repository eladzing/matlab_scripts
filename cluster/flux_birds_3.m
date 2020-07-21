%%plot t-mdot - M
%% predefine halopath smallbox/bigbox  rmfac rxfac tn tx fln flx
%% clustername pflag

tn=-1.7;tx=1.0;
fln=-0.08;flx=0.08;
smallbox=1;
bigbox=8;

rolim=[fln flx];
tlim=[tn tx];
    
    len=[200 200];

if readFlag  
    [full_ff, full_ro, full_ts, ~, rp, mvir, rvir, tvir, ~]=catspheres(smallbox,bigbox,'vcm','Rvcm','hubbleflag');
end

%% inner zone
if ~skip
    rfac=[0.01 0.2];
    
    ind=find( rp>rfac(1).*rvir & rp<rfac(2).*rvir);
    
    innerr=ind(1)-1;
    if innerr<2
        innerr=2;
    end
    outerr=ind(end)+1;
    if outerr>length(rp)
        outerr=length(rp);
    end
    clear ind
    
    % extract relavent info
    xx=full_ts(innerr:outerr,:,:);
    val=full_ro(innerr:outerr,:,:);
    yy=full_ff(innerr:outerr,:,:);
    
    [Mg_prof, ~, ~]=read_Mass_Profiles(rp(innerr-1:outerr));
    
    Mg_prof=diff(Mg_prof);%   zeros(size(Mg1));
    
    w=sum(sum(val,2),3);
    
    for i=1:length(w);
        val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
    end;
    clear w; clear Mg_prof;
    
    %% calculate 2d-histogram
    
    s=size(yy);
    omnorm=4*pi/(s(2).*s(3)); %solid angle normalization
    [mnorm,~,~]=read_Mass_Profiles(rvir); %mass norm
    %mnorm=read_MGAS_Profile(HALO_PATH, rvir); %mass norm
    tnorm=tvir;ronorm=mnorm;
    finorm=omnorm.*1e-9;
    
    % reshape
    xx=reshape(xx,numel(xx),1);
    yy=reshape(yy,numel(yy),1);
    val=reshape(val,numel(val),1);
    
    lx=log10(xx./tnorm);clear xx;
    ly=yy./ronorm./finorm;clear yy;
    
    
    % prescribed limits
        
    [troIn, binsize, xxlimIn,yylimIn]= histogram2d(lx,ly,val,...
        'xlim',tlim,'ylim',rolim,'len',len);
    %[tro bx by]=basic_bird(lx,ly,val,wt,tlim,rolim,len);
    
    binormIn=prod(binsize);
    
    cnormIn=sum(sum(squeeze(troIn(:,:,1))));
    
    %% outer
    rfac=[0.2 1.2];
    
    ind=find(rp>rfac(1).*rvir & rp<rfac(2).*rvir);
    
    innerr=ind(1)-1;
    if innerr<2
        innerr=2;
    end
    outerr=ind(end)+1;
    if outerr>length(rp)
        outerr=length(rp);
    end
    clear ind
    
    % extract relavent info
    xx=full_ts(innerr:outerr,:,:);
    val=full_ro(innerr:outerr,:,:);
    yy=full_ff(innerr:outerr,:,:);
    
    [Mg_prof, ~, ~]=read_Mass_Profiles(rp(innerr-1:outerr));
    
    Mg_prof=diff(Mg_prof);%   zeros(size(Mg1));
    
    w=sum(sum(val,2),3);
    
    for i=1:length(w);
        val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
    end;
    clear w; clear Mg_prof;
    
    %% calculate 2d-histogram
    
    s=size(yy);
    omnorm=4*pi/(s(2).*s(3)); %solid angle normalization
    [mnorm,~,~]=read_Mass_Profiles(rvir); %mass norm
    %mnorm=read_MGAS_Profile(HALO_PATH, rvir); %mass norm
    tnorm=tvir;ronorm=mnorm;
    finorm=omnorm.*1e-9;
    
    % reshape
    xx=reshape(xx,numel(xx),1);
    yy=reshape(yy,numel(yy),1);
    val=reshape(val,numel(val),1);
    
    lx=log10(xx./tnorm);clear xx;
    ly=yy./ronorm./finorm;clear yy;
    
    
    % prescribed limits
  
    
    [troOut, binsize, xxlimOut,yylimOut]= histogram2d(lx,ly,val,...
        'xlim',tlim,'ylim',rolim,'len',len);
    %[tro bx by]=basic_bird(lx,ly,val,wt,tlim,rolim,len);
    
    binormOut=prod(binsize);
    
    cnormOut=sum(sum(squeeze(troOut(:,:,1))));
    
end


%% plot bird


figure;
 %set(gcf,'Papersize',[25 15])
subplot(1,12,1:5)
imagesc(xxlimIn,yylimIn,log10(troIn(:,:,1)./cnormIn));
axis square
set(gca,'Ydir','normal','Fontsize',12);
xlabelmine('$\log\left(T/T_{vir}\right)$');
ylabelmine('$\dot{M} / M_{gas}(R_\mathrm{vir}) / d\Omega\,[\mathrm{Gyr^{-1}}]$');
grid;
%bar=colorbar;
%set(get(bar,'Title'),'String',sprintf('log(M/M_{zone})'))
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
caxis([-5 -3])

subplot(1,12,6:12)
imagesc(xxlimOut,yylimOut,log10(troOut(:,:,1)./cnormOut));
axis square
set(gca,'Ydir','normal','Fontsize',12,'YTickLabel',[]);
xlabelmine('$\log\left(T/T_{vir}\right)$');
grid;

%load('MyColormaps','avijet_bird')
%set(gcf,'Colormap',avijet_bird)

bar=colorbar;
set(get(bar,'Title'),'String','$\log(M/M_{zone})$','Interpreter','latex','Fontsize',12)
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
caxis([-5 -3])

%title(sprintf( '%s Mass Histogram (%s<r/R_{vir}<%s)',CLST,num2str(rmfac,2),num2str(rxfac,2)),'Fontsize',12);

