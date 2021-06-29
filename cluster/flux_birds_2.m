%%plot t-mdot - M
%% predefine halopath smallbox/bigbox  rmfac rxfac tn tx fln flx
%% clustername pflag

if readFlag
 [full_ff, full_ro, full_ts, ~, rp, mvir, rvir, tvir, ~]=catspheres(smallbox,bigbox,'vcm','Rvcm','hubbleflag');   
end

% choose shell of matter to examine 
if ~exist('rfac','var')
    rfac=[0 10];
end
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


%clear full_ts full_ro full_ff;

%global HALO_PATH
%Mg1=read_MGAS_Profile(HALO_PATH, rp(innerr:outerr));
[Mg_prof, ~, ~]=read_Mass_Profiles(rp(innerr-1:outerr));
   
Mg_prof=diff(Mg_prof);%   zeros(size(Mg1));
% if innerr>1
%     [Mg2,~,~]=read_Mass_Profiles(rp(innerr-1:outerr-1));
%     %Mg2=read_MGAS_Profile(HALO_PATH, rp(innerr-1:outerr-1));
%     Mg_prof=Mg1-Mg2;
% else
%     [Mg2,~,~]=read_Mass_Profiles(rp(innerr-1:outerr-1));
%     %Mg2=read_MGAS_Profile(HALO_PATH, rp(innerr:outerr-1));
%     Mg_Prof(1)=Mg1(1); Mg_prof(2:length(Mg_prof))=Mg1(2:length(Mg1))-Mg2;
% end
% clear Mg1 Mg2

w=sum(sum(val,2),3);

for i=1:length(w);
val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
end;
clear w; clear Mg_prof;

%weight=1; %yyal<0;
%wt=ones(size(xx));
%wt=wt.*weight;

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
%%val=val%%./mnorm;

% prescribed limits
rolim=[fln flx];
tlim=[tn tx];

len=[500 500];

[tro, binsize, xxlim,yylim]= histogram2d(lx,ly,val,...
    'xlim',tlim,'ylim',rolim,'len',len);
%[tro bx by]=basic_bird(lx,ly,val,wt,tlim,rolim,len);

binorm=prod(binsize);

%% plot bird

cnorm=squeeze(tro(:,:,1));
figure; imagesc(tlim,rolim,log10(tro(:,:,1)./sum(cnorm(:))));
bar=colorbar;
set(get(bar,'Title'),'String',sprintf('log(M/M_{zone});pixel=%s',num2str(binorm,'%1.2g')),'Fontsize',12)
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
caxis([-6.2 -2])
%set(gca,'Xdir','normal');
set(gca,'Ydir','normal');
xlabelmine('$\log\left(T/T_{vir}\right)$');
ylabelmine('$\dot{M} / Mgas_{vir} / d\Omega\,[\mathrm{Gyr^{-1}}]$');
grid;
%title(sprintf( '%s Mass Histogram (%s<r/R_{vir}<%s)',CLST,num2str(rmfac,2),num2str(rxfac,2)),'Fontsize',12);
set(gca,'fontsize',12)

% %% print to file
% if strcmp(pflag,'print')
%   
%   %result_dir='/home/titan3/eladzing/cold_flows/printout';
%   saveas(gcf,sprintf('%s/%s_tfluxbird_%s.png',result_dir,CLST,printag));
% end

% %% plot hist
% figure;
% mdhist=sum(tro(:,:,1),2); %create 
% 
% dt=binsize(2);
% tax=fln:dt:flx-dt;
% 
% plot(tax,mdhist./dt); grid; xlim([fln flx]);
% set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
% xlabelmine('$\dot{M}/Mgas_{vir}/d\Omega \,[\mathrm{Gyr^{-1}}]$');
% ylabelmine('$\frac{d(M/Mgas)} {d\log(\dot{M})}$')
% %title(sprintf( '%s Flux Mass Histogram (%s<r/R_{vir}<%s)',CLST,num2str(rmfac,2),num2str(rxfac,2)),'fontsize',12);
% set(gca,'fontsize',12)
% %%legend('Inflow','Outflow','Location','NorthWest');
% 
% % if strcmp(pflag,'print')
% %   
% %   %result_dir='/home/titan3/eladzing/cold_flows/printout';
% %   saveas(gcf,sprintf('%s/%s_fluxhist_%s.png',result_dir,CLST,printag));
% % end
% 
% 
% 
