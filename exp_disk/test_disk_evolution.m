rstrip=1.5;
fg=0.25;
beta=0.5;

ls=9:0.1:12;
ms=10.^ls;

rd=1:0.1:5;

[msMat rdMat]=meshgrid(ms,rd);

sigMat=(0.5.*msMat)./(pi.*(1.67834.*rdMat).^2);

sigNew=zeros(size(sigMat));

for i=1:length(ms)
    for j=1:length(rd)
        disc=mk_disk_struct(ms(i),rd(j),fg,beta);
        [sigma50,~]=adjust_exp_disk_sf_old(disc,rstrip.*disc.rd);
        sigNew(j,i)=sigma50;
    end
end
fac=sigNew./sigMat;
[sigma501,~]=adjust_exp_disk_sf(rstrip,fg,beta);




% imagesc(ls,rd,log10(sigMat))
% set(gca,'Ydir','normal','Fontsize',12);
% bar=colorbar;








% %% testing for beta
% hf(1)=figure;
% hf(2)=figure;
% rs0=0:0.01:20;
% h1=[];h2=[];
% c={'r' 'b' 'g'};
% b=[0.5 1 2];
% for j=1:length(b);
% disc=mk_disk_struct(1e10,5,0.1,b(j));
% s1=zeros(size(rs0));
% s2=s1;
% %r50=s1;
% %reff=effective_rad(0.5);
% rs=rs0.*disc.rd;
% for i=1:length(rs);
% [s1(i),s2(i),~]=adjust_exp_disk_sf(disc,rs(i));
% end
% figure(hf(1));
% h1(end+1)=plot(rs./disc.rd,s1./disc.sigma_eff,c{j},'LineWidth',2);
% hold on
%
% figure(hf(2));
% h2(end+1)=plot(rs./disc.rd,s2./disc.sigma_5090,c{j},'LineWidth',2);
% hold on
% end
% figure(hf(1));
% hl=legend(h1,'$\beta=0.5$','$\beta=1$','$\beta=2$');
% set(hl,'Interpreter','latex','Location','SouthEast')
% xlabelmine('$r_{strip}/R_d$')
% ylabelmine('$\Sigma_{\mathrm{eff}}/\Sigma_{\mathrm{eff}}^0$')
% grid
%
% figure(hf(2));
% hl=legend(h2,'$\beta=0.5$','$\beta=1$','$\beta=2$');
% set(hl,'Interpreter','latex','Location','SouthWest')
% xlabelmine('$r_{strip}/R_d$')
% ylabelmine('$\Sigma_{5090}/\Sigma_{5090}^0$')
% grid

% %% testing for fg & Beta
% hf(1)=figure;
% hf(2)=figure;
% rs0=0:0.01:20;
% h1=[];h2=[];
% c={'r' 'b' 'g'};
% lt={'--','-',':'};
% b=[0.5 1 2];
% f=[0.1 0.25 0.5];
% for j=1:length(b);
%     for k=1:length(f)
%         disc=mk_disk_struct(1e10,5,f(k),b(j));
%         s1=zeros(size(rs0));
%         s2=s1;
%         %r50=s1;
%         %reff=effective_rad(0.5);
%         rs=rs0.*disc.rd;
%         for i=1:length(rs);
%             [s1(i),s2(i),~]=adjust_exp_disk_sf(disc,rs(i));
%         end
%         ll=strcat(c{j},lt{k});
%         figure(hf(1));
%         h1(end+1)=plot(rs./disc.rd,s1./disc.sigma_eff,ll,'LineWidth',2,...
%         'DisplayName',sprintf('$f_g=%s,\\beta=%s$',num2str(f(k)),num2str(b(j))));
%         hold on
%
%         figure(hf(2));
%         h2(end+1)=plot(rs./disc.rd,s2./disc.sigma_5090,ll,'LineWidth',2,...
%          'DisplayName',sprintf('$f_g=%s,\\beta=%s$',num2str(f(k)),num2str(b(j))));
%         hold on
%     end
% end
% figure(hf(1));
% hl=legend(h1);
% set(hl,'Interpreter','latex','Location','NorthEastOutside','fontsize',12)
% xlabelmine('$r_{strip}/R_d$')
% ylabelmine('$\Sigma_{\mathrm{eff}}/\Sigma_{\mathrm{eff}}^0$')
% grid
% set(gca,'fontsize',12)
%
% figure(hf(2));
% hl=legend(h2);
% set(hl,'Interpreter','latex','Location','NorthEastOutside','fontsize',12)
% xlabelmine('$r_{strip}/R_d$')
% ylabelmine('$\Sigma_{5090}/\Sigma_{5090}^0$')
% grid
% set(gca,'fontsize',12)