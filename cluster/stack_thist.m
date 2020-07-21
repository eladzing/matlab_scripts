
stv1=zeros(slen-1,size(stack{1,3},1));
stv2=stv1;
for i=1:slen
stack{i,1}
stv1(i,:)=stack{i,3};
stv2(i,:)=stack{i,4};
end
vmin=mean(stv1,1);
vmout=mean(stv2,1);

figure;plot(tax,[vmin' vmout']'); grid; xlim([mnt mxt]);
xlabel('log T/T_{vir}');
ylabel('d(M/Mgas)/dlog(T/T_{vir})');
legend('Inflow', 'Outflow','Location','NorthWest');
title(sprintf( '%s Infolw/Outflow Temperature Histogram (r/R_{vir}<0.02)','Stacked'));
saveas(gcf,sprintf('%s/%s_thist_core.png',result_dir,'stack'));
