new_env(CL);
ss=S(box);
ro=RHOG(box);

lss=log10(ss);lro=log10(ro);
%%%mnro=min(lro(:));mxro=max(lro(:));
%%mnss=min(lss(:));mxss=max(lss(:));
mnro=7;
mxro=17;
mnss=-8;
mxss=3;


[sro bx by]=basic_bird(lss,lro,ro,ones(size(ro)),[mnss mxss],[mnro mxro],[0 0]);


cnorm=squeeze(sro(:,:,1));
figure; imagesc([mnss mxss],[mnro mxro],log10(sro(:,:,1)./sum(cnorm(:))));
load('MyColormaps','avijet_bird')
set(gcf,'Colormap',avijet_bird)
set(gca,'Ydir','normal');
xlabel('log S','Fontsize',14);
ylabel('log \rho','Fontsize',14);
set(gca,'fontsize',12)
bar=colorbar;

set(get(bar,'Title'),'String','log(M/M_{gas})','Fontsize',12)
result_dir='/home/titan3/eladzing/cold_flows/printout';
title(sprintf( '%s Mass Histogram (box=%d Mpc/h',CL,box),'Fontsize',14);
saveas(gcf,sprintf('%s/%s_sro_bird_%d.png',result_dir,CL,box));