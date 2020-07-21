%%plot t-mdot - M

[full_ff full_ro full_ts rp m200 r200 t200]=catspheres(halopath,smallbox,bigbox);   
    
ind=find( rp>rmfac.*r200 & rp<rxfac.*r200);
innerr=ind(1)-1;
outerr=ind(length(ind))+1;

%innerr=; outerr=length(rp);

xx=full_ts(innerr:outerr,:,:);
val=full_ro(innerr:outerr,:,:);
yy=full_ff(innerr:outerr,:,:);
clear full_ts full_ro full_ff;

Mg_prof=read_MGAS_Profile(halopath, rp(innerr:outerr));

w=sum(sum(val,2),3);

for i=1:length(w);
val(i,:,:)=val(i,:,:).*Mg_prof(i)./w(i);
end;
clear w; clear Mg_prof;

weight=1; %yy<0;
wt=ones(size(xx));
wt=wt.*weight;
%% calculate 2d-histogram
s=size(yy);
tnorm=log10(t200);ronorm=4*pi/(s(2).*s(3));

ltm=xx(:);clear xx;
lro=yy(:);clear yy;
wt=wt(:);
val=val(:);

ltm=log10(ltm);
%lro=log10(abs(lro));

% prescribed limits
romin=fln;romax=flx;
indd=find((lro>romin)&(lro<romax));
ltm=ltm(indd);
lro=lro(indd);
%lro=lro.*mask;
%wt=wt.*mask;
%val=val.*mask;

len=200;
mnro=romin;% floor(min(lro));%;  min(min(lro))));
%mnro=floor(min(lro));%;  min(min(lro))));
mxro=romax;%ceil(max(lro));%-mnro;%  max(max(lro))));
%mxro=ceil(max(lro));%-mnro;%  max(max(lro))));
mnt=floor(min(ltm));%min(min(ltm))));
mxt=ceil(max(ltm));%-mnt;  %max(max(ltm))));

tro=zeros(len,len,2);

for i=1:size(ltm,1)   %.*size(ltm,2).*size(ltm,3))
    indx=ceil(((ltm(i)-mnt)./(mxt-mnt)).*len);
    indy=ceil(((lro(i)-mnro)./(mxro-mnro)).*len);  %len+1-ceil((lro(i)./mxro).*len);  
    if(indx>=1 && indx<=len) && (indy>=1 && indy<=len)
        tro(indy,indx,1)=tro(indy,indx,1)+val(i).*wt(i);
        tro(indy,indx,2)=tro(indy,indx,2)+wt(i);
    end
end

%% plot bird
%ind=find(tro(:,:,2)==0);
%tro(ind)=1e-50;

figure; imagesc([mnt-tnorm mxt-tnorm],[mnro/ronorm mxro/ronorm],log10(tro(:,:,1)));
bar=colorbar;
set(get(bar,'Title'),'String','M [M_{sun}]')
load('MyColormaps','birdcmap')
set(gcf,'Colormap',birdcmap)
%set(gca,'Xdir','normal');
set(gca,'Ydir','normal');
xlabel('log T/T_{vir}');
ylabel('Mdot/solid angle [M_{sun}/yr/rad^2]');
title(sprintf( '%s Mass Histogram (0.02<r/R_{vir}<0.2)',clustername));
%set(gca,'XTick',[(2-mnt)*len/mxt,(3-mnt)*len/mxt,(4-mnt)*len/mxt,(5-mnt)*len/mxt,(6-mnt)*len/mxt,(7-mnt)*len/mxt]);
%set(gca,'YTick',[(-30-mnro)*len/mxro,(-28-mnro)*len/mxro,(-26-mnro)*len/mxro,(-24-mnro)*len/mxro,(-22-mnro)*len/mxro,(-20-mnro)*len/mxro,(-18-mnro)*len/mxro]);
%set(gca,'XTickLabel',[2,3,4,5,6,7]);
%set(gca,'YTickLabel',[-30,-28,-26,-24,-22,-20,-18]);

%% print to file
if strcmp(pflag,'print')
  
  result_dir='/home/titan3/eladzing/cold_flows/printout';
  saveas(gcf,sprintf('%s/%s_tfluxbird_core.png',result_dir,clustername));
end
  
