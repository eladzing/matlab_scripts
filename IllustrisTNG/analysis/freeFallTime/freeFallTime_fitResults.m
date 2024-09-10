global illUnits 
global cosmoStruct
global simDisplayName

%load('/home/zinger/workProjects/matlab_scripts/IllustrisTNG/matFiles/freeFallTime_profiles_snp99_TNG100.mat')

mv=tffProfile.NFW.mv(mask);
galMass=subs.SubhaloMassInRadType(illustris.partTypeNum('stars')+1,mask).*illUnits.massUnit;

aCof=tffProfile.polyfit.a(mask);
bCof=tffProfile.polyfit.b(mask);
cCof=tffProfile.polyfit.c(mask);

cof=cat(1,aCof,bCof,cCof);

lr=(-2:0.02:log10(2))';

lrr=cat(2,lr.^2,lr,ones(size(lr)));

profs=lrr*cof;

% myFigure;
% plot(lr,profs)


%% plot profiles in stellar mass bins 
fnam='tffProfile_mstRange_%s_%s_snp%s_%s';

lm=log10(galMass);
ll=cat(1,lr,flipud(lr));
massB=[9:0.5:11.5 Inf];
cc=brewermap(8,'Set1');


for i=1:length(massB)-1
    
    binMask=lm>=massB(i) & lm<massB(i+1);
    
    pr=profs(:,binMask);
    mpr(:,i)=median(pr,2);
    quant(:,:,i)=quantile(pr,[0.1 0.9],2);
    
        
    myFigure;
    plot(lr,pr,'color',[0.75 0.75 0.75])
    hold on
    
    bb=cat(1,quant(:,1,i),flipud(quant(:,2,i)));
    patch(ll,bb,cc(1,:),'faceAlpha',0.25,'edgecolor','none')
      
    plot(lr,mpr(:,i),'color',cc(2,:),'linewidth',2);
    
    set(gca,'fontsize',14);
    grid   
    titlemine(sprintf('halo mass range:%s - %s',num2str(massB(i)),num2str(massB(i+1))));
    xlabelmine('$\log\, r/R_\mathrm{200,c}$');
    ylabelmine('$\log t_\mathrm{ff} \, [\mathrm{Gyr}]$');
    
   % printout_fig(gcf,...
   %     sprintf(fnam,num2str(massB(i)),num2str(massB(i+1)),num2str(snap),simDisplayName))
end

%%


h=[];
myFigure;

kCol=[1:5 7];
for i=1:length(massB)-1
    nam=[num2str(massB(i)) '-' num2str(massB(i+1))];
    
    bb=cat(1,quant(:,1,i),flipud(quant(:,2,i)));
    patch(ll,bb,cc(kCol(i),:),'faceAlpha',0.25,'edgecolor','none')
    
    h(i)= plot(lr,mpr(:,i),'color',cc(kCol(i),:),'DisplayName',nam,'linewidth',2);
    if i==1
        hold on
    end
    
end

hl=legend(h); %'9-9.5','9.5-10','10-10.5','10.5-11','11-11.5','11.5<')
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
grid
set(gca,'fontsize',14);
xlabelmine('$\log\, r/R_\mathrm{200,c}$');
ylabelmine('$\log t_\mathrm{ff} \, [\mathrm{Gyr}]$');  
titlemine('$t_\mathrm{ff}$ profiles in stellar mass bins');
 printout_fig(gcf,...
        sprintf('tffProfile_mst_%s_%s_snp%s_%s',num2str(snap),simDisplayName));


%% plot profiles in halo mass bins 
fnam='tffProfile_mhalRange_%s_%s_snp%s_%s';

lm=log10(mv);
ll=cat(1,lr,flipud(lr));
cc=brewermap(8,'Set1');
massB=[10:1:14 Inf];


for i=1:length(massB)-1
    
    binMask=lm>=massB(i) & lm<massB(i+1);
    
    pr=profs(:,binMask);
    mpr(:,i)=median(pr,2);
    quant(:,:,i)=quantile(pr,[0.1 0.9],2);
    
        
    myFigure;
    plot(lr,pr,'color',[0.75 0.75 0.75])
    hold on
    
    bb=cat(1,quant(:,1,i),flipud(quant(:,2,i)));
    patch(ll,bb,cc(1,:),'faceAlpha',0.25,'edgecolor','none')
      
    plot(lr,mpr(:,i),'color',cc(2,:),'linewidth',2);
    
    %set(gca,'fontsize',14);
    grid   
    titlemine(sprintf('Halo Mass range:%s - %s',num2str(massB(i)),num2str(massB(i+1))));
    xlabelmine('$\log\, r/R_\mathrm{200,c}$');
    ylabelmine('$\log t_\mathrm{ff} \, [\mathrm{Gyr}]$');
    myAxis; 
 %   printout_fig(gcf,...
 %       sprintf(fnam,num2str(massB(i)),num2str(massB(i+1)),num2str(snap),simDisplayName))

    
    
end


h=[];
myFigure;
kCol=1:5;
for i=1:length(massB)-1
    nam=[num2str(massB(i)) '-' num2str(massB(i+1))];
    
    bb=cat(1,quant(:,1,i),flipud(quant(:,2,i)));
    patch(ll,bb,cc(kCol(i),:),'faceAlpha',0.25,'edgecolor','none')
    
    h(i)= plot(lr,mpr(:,i),'color',cc(kCol(i),:),'DisplayName',nam,'linewidth',2);
    if i==1
        hold on
    end
    
end

hl=legend(h); %'9-9.5','9.5-10','10-10.5','10.5-11','11-11.5','11.5<')
set(hl,'Interpreter','latex','fontsize',14,'location','SouthEast')
grid
set(gca,'fontsize',14);
xlabelmine('$\log\, r/R_\mathrm{200,c}$');
ylabelmine('$\log t_\mathrm{ff} \, [\mathrm{Gyr}]$');  
titlemine('$t_\mathrm{ff}$ profiles in halo mass bins');
 %printout_fig(gcf,...
 %       sprintf('tffProfile_mhal_%s_%s_snp%s_%s',num2str(snap),simDisplayName));



%% plot residulas of fit 
for i=1:6
ptsP(i) = mk_meanMedian_bin(log10(mv),tffProfile.polyfit.residuals(i,mask),'nb',20);
end


myFigure;
cc=brewermap(6,'Set1');
for k=1:6
    nam=num2str(k);
    h(k)=plot(ptsP(k).xMedian,ptsP(k).yMedian,'color',cc(k,:),'Displayname',nam);
    if k==1
        hold on
    end
    
    plot(ptsP(k).xMedian,ptsP(k).yQuarts(1,:),':','color',cc(k,:));
    plot(ptsP(k).xMedian,ptsP(k).yQuarts(4,:),':','color',cc(k,:));
    
end
hl=legend(h);
xlabelmine('$\log\,M_{200}$');
ylabelmine('residuals %');
titlemine('polyfit');

for i=1:6
ptsN(i) = mk_meanMedian_bin(log10(mv),tffProfile.NFW.residuals(i,mask),'nb',20);
end


myFigure;
cc=brewermap(6,'Set1');
for k=1:6
    nam=num2str(k);
    h(k)=plot(ptsN(k).xMedian,ptsN(k).yMedian,'color',cc(k,:),'Displayname',nam);
    if k==1
        hold on
    end
    
    plot(ptsN(k).xMedian,ptsN(k).yQuarts(1,:),'--','color',cc(k,:));
    plot(ptsN(k).xMedian,ptsN(k).yQuarts(4,:),'--','color',cc(k,:));
    
end
hl=legend(h);
xlabelmine('$\log\,M_{200}$');
ylabelmine('residuals %');
titlemine('NFW')
