%%plot t-mdot - M
%% predefine halopath smallbox/bigbox  rmfac rxfac tmfac
%% clustername pflag
[full_ff full_ro full_ts rp m200 r200 t200]=catspheres(halopath,smallbox,bigbox);   
clear full_ff

%tmfac=0.333333;
%rxfac=0.5;
%rmfac=0.0;
ind=find( rp>rmfac.*r200 & rp<rxfac.*r200);
innerr=ind(1)-1; 

if innerr<1 
    innerr=1;
end
outerr=ind(length(ind));
if outerr>=length(rp)
    outerr=length(rp)-1;
end

tm=full_ts(innerr:outerr,:,:);
val=full_ro(innerr:outerr,:,:);
%test=full_ff(innerr:outerr,:,:);
clear full_ts full_ro ; 
%tnorm=t200;

%mnt=floor(log10(min(tm(:))./tnorm));%min(min(ltm))));
%mxt=ceil(log10(max(tm(:))./tnorm));

%% change density values to mass values. 

%rp2=rp1;
if innerr>1
  rp1= 0.5.*(rp(innerr:outerr)+rp(innerr-1:outerr-1));    %% rp2 is the radius halfway between bins
else
  rp1(1)=0.5.*rp(1);
  rp1(2:outerr)= 0.5.*(rp(2:outerr)+rp(1:outerr-1));   %% rp2 is the radius halfway between bins
end

rp2= 0.5.*(rp(innerr+1:outerr+1)+rp(innerr:outerr));

Mg1= read_MGAS_Profile(halopath,rp1);
Mg2= read_MGAS_Profile(halopath,rp2);
Mg_prof=Mg2-Mg1;
clear rp1 rp2 Mg1 Mg2;

%mnorm=sum(Mg_prof);  %(length(Mg_prof))-Mg_prof(1);

w=sum(sum(val,2),3);


for i=1:length(w);
val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
end;
clear w; clear Mg_prof;

%% sum up in log bins


mask=log10(tm/t200)<(tmfac);
lfac=10^(tmfac);
mval=val.*mask;

tfprof1=transpose(cumsum(sum(sum(mval,2),3)));
mnorm1=transpose(cumsum(sum(sum(val,2),3)));

lthick=0.05;
lbin=-3:lthick:log10(rp(outerr)./r200);
inr=10.^lbin;
rbin=10.^(lbin(1:length(lbin)-1)+0.5.*lthick);
tfprof2=[];
mnorm2=[];

for j=1:length(inr)-1
    
    ind=find((rp(innerr:outerr)>inr(j).*r200) & (rp(innerr:outerr)<inr(j+1).*r200));
    t1=mval(ind,:,:);
    m1=val(ind,:,:);
    tfprof2(end+1)=sum(t1(:));
    mnorm2(end+1)=sum(m1(:));
    
    if i>outerr; break; end
end


clear tm mask val ;
 
%% plot 

if strcmp(plotflag,'plot')
    figure;subplot(2,1,1);
    loglog(rp(innerr:outerr)./r200,(tfprof1./mnorm1)); grid;
    %set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
    xlabel('r/R_{vir}');
    ylabel(sprintf('M(<r,T<%sT_{vir})/M(<r)',num2str(lfac,2)));
    title(sprintf( '%s Cumulative and Total Cold Gas fraction (T<%sT_{vir}) Profile',clustername,num2str(lfac,2)));
    %legend('Inflow', 'Outflow','Location','NorthWest');

    subplot(2,1,2);
    semilogx(rbin,tfprof2./mnorm2,'.-');grid;
    %set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
    xlabel('r/R_{vir}');
    ylabel(sprintf('M(T<%sT_{vir},r)/M(r)',num2str(lfac,2)));

end





%bar=colorbar;
%set(get(bar,'Title'),'String','log M/Mgas_{vir}')
%load('MyColormaps','birdcmap')
%set(gcf,'Colormap',birdcmap)
%set(gca,'Xdir','normal');
%set(gca,'Ydir','normal');

%set(gca,'XTick',[(2-mnt)*len/mxt,(3-mnt)*len/mxt,(4-mnt)*len/mxt,(5-mnt)*len/mxt,(6-mnt)*len/mxt,(7-mnt)*len/mxt]);
%set(gca,'YTick',[(-30-mnro)*len/mxro,(-28-mnro)*len/mxro,(-26-mnro)*len/mxro,(-24-mnro)*len/mxro,(-22-mnro)*len/mxro,(-20-mnro)*len/mxro,(-18-mnro)*len/mxro]);
%set(gca,'XTickLabel',[2,3,4,5,6,7]);
%set(gca,'YTickLabel',[-30,-28,-26,-24,-22,-20,-18]);

%% print to file
if strcmp(pflag,'print')
    
  saveas(gcf,sprintf('%s/%s_tfrac.png',result_dir,clustername));
end
  
