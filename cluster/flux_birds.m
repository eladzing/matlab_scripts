%%plot t-ro -mdot

[full_ff full_ro full_ts rp m200 r200 t200]=catspheres(halopath,smallbox,bigbox);   


ind=find(rp>0.01.*RVIR & rp<0.1*RVIR);
innerr=ind(1)-1;
outerr=ind(length(ind))+1;
%innerr=1; outerr=length(rp);

xx=full_ts(innerr:outerr,:,:);
yy=full_ro(innerr:outerr,:,:);
val=ones(size(xx));full_ff(innerr:outerr,:,:);

weight=1;

%% calculate 2d-histogram
ltm=log10(xx);
lro=log10(yy);

wt=ones(size(ltm));
wt=wt.*weight;

ltm=ltm(:);
lro=lro(:);
wt=wt(:);
val=val(:);

len=200;
mnro=floor(min(lro));%  min(min(lro))));
mxro=ceil(max(lro));%  max(max(lro))));
mnt=floor(min(ltm)); %min(min(ltm))));
mxt=ceil(max(ltm));  %max(max(ltm))));

tro=zeros(len,len,2);
%ltm=ltm-mnt;
%lro=lro-mnro;
for i=1:size(ltm,1)   %.*size(ltm,2).*size(ltm,3))
    indx=ceil(((ltm(i)-mnt)./mxt).*len);
    indy=ceil(((lro(i)-mnro)./mxro).*len);  %len+1-ceil((lro(i)./mxro).*len);  
    tro(indy,indx,1)=tro(indy,indx,1)+val(i).*wt(i);
    tro(indy,indx,2)=tro(indy,indx,2)+wt(i);
end

%% plot bird
%ind=find(tro(:,:,2)==0);
%tro(ind)=1e-50;

figure; imagesc(log10(tro(:,:,1)));
colorbar;
load('MyColormaps','birdcmap')
set(gcf,'Colormap',birdcmap)
set(gca,'Xdir','normal');
set(gca,'Ydir','normal');
xlabel('log T');
ylabel('log \rho');
title(sprintf('%s Flux Histogram on T-\rho Plane',clustername);
%set(gca,'XTick',[(2-mnt)*len/mxt,(3-mnt)*len/mxt,(4-mnt)*len/mxt,(5-mnt)*len/mxt,(6-mnt)*len/mxt,(7-mnt)*len/mxt]);
%set(gca,'YTick',[(-30-mnro)*len/mxro,(-28-mnro)*len/mxro,(-26-mnro)*len/mxro,(-24-mnro)*len/mxro,(-22-mnro)*len/mxro,(-20-mnro)*len/mxro,(-18-mnro)*len/mxro]);
%set(gca,'XTickLabel',[2,3,4,5,6,7]);
%set(gca,'YTickLabel',[-30,-28,-26,-24,-22,-20,-18]);

%% print to file
if pflag=='print'
  
  result_dir='/home/titan3/eladzing/cold_flows/printout';
  saveas(gcf,sprintf('%s/%s_rotee_fluxbird.png',result_dir,clustername));
end
