

%%set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
global HALO_PATH
%%load(sprintf('%s/virial%d', HALO_PATH, 8));

[full_ff full_ro full_ts full_vr rp m200 r200 t200]=catspheres(HALO_PATH,smallbox,bigbox);
clear full_ff;
ind=find( rp>rmfac.*r200 & rp<rxfac.*r200);

innerr=ind(1)-1; 
if innerr<1 
    innerr=1;
end
outerr=ind(length(ind))+1;
if outerr>length(rp)
    outerr=length(rp);
end

fts=full_ts(innerr:outerr,:,:);
fro=full_ro(innerr:outerr,:,:);
val=fro;
fvr=full_vr(innerr:outerr,:,:);
rps=rp(innerr:outerr);
clear full_ts full_ro full_vr;

% create mass sphere
Mg1=read_MGAS_Profile(HALO_PATH, rp(innerr:outerr));
Mg_prof=zeros(size(Mg1));
if innerr>1
    Mg2=read_MGAS_Profile(HALO_PATH, rp(innerr-1:outerr-1));
    Mg_prof=Mg1-Mg2;
else
    Mg2=read_MGAS_Profile(HALO_PATH, rp(innerr:outerr-1));
    Mg_Prof(1)=Mg1(1); Mg_prof(2:length(Mg_prof))=Mg1(2:length(Mg1))-Mg2;
end
clear Mg1 Mg2

w=sum(sum(val,2),3);

for i=1:length(w);
val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
end;
clear w; clear Mg_prof;

Ggrav=4.43e-15;
%zmet=1.;

%[tlam zlam lambda]=read_lambda; %%read interpolation data for cooling function (in log10) 

zmet=1*ones(size(fts));
 
%lamb=interp2(zlam,tlam,lambda,zmet,log10(fts),'spline');
%lamb22=10.^(lamb-(-22));


%clear lamb zlam tlam lambda; 

lmb22=0.6.*(zmet/0.3).^0.7.*(fts./1e6).^(-1)+0.02.*(fts./1e6).^0.5;

tcool=2.61.*(fro.*6.77e-13).^(-1).*(fts./1e6).*lmb22.^-1; %in Gyr
%%tdyn2=(4.*pi.*Ggrav.*fro).^(-1/2);
cson=sqrt(1.52e8.*fts); %speed of sound in cm/sec


fullr=ones(size(fvr));

for i=length(rps)
    fullr(i,:,:)=rps(i);
end

%tcon=fullr./abs(fvr).*(979.365); % in Gyr
tdcs=fullr./cson.*3.085e24./(3.15e16);


clear fullr fvr fts fro rps rp; 

%tlim=[-4 4];
tlimx=[-4 4];tlimy=[-5 5];
len=[200 200];

[tro bx by]=basic_bird(log10(tcool),log10(tdcs),val,ones(size(val)),tlimx,tlimy,len);
binorm=bx.*by;
if strcmp(plotflag,'plot')
    cnorm=squeeze(tro(:,:,1));
    figure; imagesc(tlimx,tlimy,log10(tro(:,:,1)./sum(cnorm(:))));
    load('MyColormaps','birdmaxcmap')
    set(gcf,'Colormap',birdmaxcmap)
    bar=colorbar;
    set(get(bar,'Title'),'String',sprintf('log(M/M_{zone});pixel=%s',num2str(binorm,'%1.2g')),'Fontsize',12)
    %caxis([-7 -2.5])
    %set(gca,'Xdir','normal');
    set(gca,'Ydir','normal');
    thub=log10([14 14]);
    line(tlimx,thub,'Color','Black');
    line(thub,tlimy,'Color','Black');
    line(tlimx,tlimx,'Color','Black');
    xlabel('log t_{cool} [Gyr]','Fontsize',12);
    ylabel('log (r/c_s) [Gyr]','Fontsize',12);grid;
    title(sprintf( '%s Mass Histogram (%s<r/R_{vir}<%s)',clustername,num2str(rmfac,2),num2str(rxfac,2)),'Fontsize',12);
    set(gca,'fontsize',12)
      
    if strcmp(pflag,'print')
        saveas(gcf,sprintf('%s/%s_timebirds2_%s.png',result_dir,clustername,printag));
    end

end


    
    
%clear all 
