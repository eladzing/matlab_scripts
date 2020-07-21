%%plot t-mdot - M
%% predefine halopath smallbox/bigbox  rmfac rxfac tn tx fln flx
%% clustername pflag
[full_ff full_ro full_ts full_vr full_ss rp m200 r200 t200 v200]=catspheres(halopath,smallbox,bigbox);   
    
ind=find( rp>rmfac.*r200 & rp<rxfac.*r200);
innerr=ind(1)-1; 
if innerr<1 
    innerr=1;
end
outerr=ind(length(ind));
if outerr>=length(rp)
    outerr=length(rp)-1;
end

vr=full_vr(innerr:outerr,:,:);
val=full_ro(innerr:outerr,:,:);
test=full_vr(innerr:outerr,:,:);
clear full_ff full_ro full_ts full_vr full_ss;
%tnorm=t200;
vnorm=v200;

%mnt=floor(log10(min(tm(:))./tnorm));%min(min(ltm))));
%mxt=ceil(log10(max(tm(:))./tnorm));

%Mg1=read_MGAS_Profile(halopath, rp(innerr:outerr));
%Mg_prof=zeros(size(Mg1));
%if innerr>1
%    Mg2=read_MGAS_Profile(halopath, rp(innerr-1:outerr-1));
%    Mg_prof=Mg1-Mg2;
%else
%    Mg2=read_MGAS_Profile(halopath, rp(innerr:outerr-1));
%    Mg_Prof(1)=Mg1(1); Mg_prof(2:length(Mg_prof))=Mg1(2:length(Mg1))-Mg2;
%end
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


mnorm=sum(Mg_prof);  %(length(Mg_prof))-Mg_prof(1);


w=sum(sum(val,2),3);

for i=1:length(w);
val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
end;
clear w; clear Mg_prof;

vr_in=abs(vr(test<0));
vr_out=abs(vr(test>0));
val_in=val(test<0);
val_out=val(test>0);
clear tm test;
%% calculate 2d-histogram

vr_in=log10(vr_in(:)./vnorm);
vr_out=log10(vr_out(:)./vnorm);
val_in=val_in(:)./mnorm;
val_out=val_out(:)./mnorm;

% prescribed limits
%tmin=tn;tmax=tx;
%indd=find((ltm>tmin)&(ltm<tmax));
%ltm=ltm(indd);
%lro=lro(indd);

%len=1e4;

mnt=vrmin;
mxt=vrmax;


vr_in=zeros(len,2);
%indx=ceil(((tm_in-mnt)./(mxt-mnt)).*len);

for i=1:size(vr_in,1)   %.*size(lvr,2).*size(lvr,3))
    indx=ceil((vr_in(i)-mnt)./(mxt-mnt).*len);
    if(indx>=1 && indx<=len)
        vr_in(indx,1)=vr_in(indx,1)+val_in(i);
        vr_in(indx,2)=vr_in(indx,2)+1;
    end
end
clear indx val_in vr_in

vr_out=zeros(len,2);
%indx=ceil(((vr_out-mnt)./(mxt-mnt)).*len)
for i=1:size(vr_out,1)   %.*size(lvr,2).*size(lvr,3))
    indx=ceil((vr_out(i)-mnt)./(mxt-mnt).*len);
    if(indx>=1 && indx<=len)
        vr_out(indx,1)=vr_out(indx,1)+val_out(i);
        vr_out(indx,2)=vr_out(indx,2)+1;
    end
end

clear indx val_out vr_out


%% plot bird

dt=(mxt-mnt)/len;
tax=mnt:dt:mxt-dt;
figure;plot(tax,[vr_in(:,1)/dt vr_out(:,1)/dt]'); grid; xlim([mnt mxt]);
set(gcf,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 1 0])
xlabel('log T/T_{vir}');
ylabel('d(M/Mgas)/dlog(T/T_{vir})');
title(sprintf( '%s Infolw/Outflow Temperature Histogram (%s<r/R_{vir}<%s)',clustername,num2str(rmfac,2),num2str(rxfac,2)));
legend('Inflow', 'Outflow','Location','NorthWest');

%% print to file
if strcmp(pflag,'print')
  
  %result_dir='/home/titan3/eladzing/cold_flows/printout';
  saveas(gcf,sprintf('%s/%s_vr_hist_%s.png',result_dir,clustername,printag));
end
  
