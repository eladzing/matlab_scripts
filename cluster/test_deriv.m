outpath='~/work/clusters/printout';

for ii=1:16
    
cname=coolprofs{ii,1};
qb=coolprofs{ii,6};
rq=qb(1,:);
q=qb(2,:);
lq=log10(q);lrq=log10(rq);

x=lrq(lrq<=-1);
y=lq(lrq<=-1);


dy1=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2));
x1=x(2:end-1);

dy2=(dy1(3:end)-dy1(1:end-2))./(x1(3:end)-x1(1:end-2));
x2=x1(2:end-1);

xl1=x(y==min(y));
xl2=x2(dy2==max(dy2));

xm=min(x);xx=max(x);    
ym=min(y);yx=max(y);    
y1m=min(dy1);y1x=max(dy1);    
y2m=min(dy2);y2x=max(dy2);    

figure;
subplot(3,1,1)

plot(x,y)
hold on
plot([xl1 xl1],[ym yx],'--k')
plot([xl2 xl2],[ym yx],'--k')
hold off
title(sprintf('%s Inner Cooling Profile Derivative Test',cname));
ylabel('log(q)');
xlim([xm xx])
subplot(3,1,2)
plot(x1,dy1)
hold on
plot([xl1 xl1],[y1m y1x],'--k')
plot([xl2 xl2],[y1m y1x],'--k')
hold off
ylabel('1^{st} deriv');
xlim([xm xx])
subplot(3,1,3)
plot(x2,dy2)
hold on
plot([xl1 xl1],[y2m y2x],'--k')
plot([xl2 xl2],[y2m y2x],'--k')
hold off
ylabel('2^{nd} deriv');
xlim([xm xx])
xlabel('log(r) [Mpc]')

saveas(gcf,sprintf('%s/%s_coolprof_testderiv.png',outpath,cname));

close gcf
end
