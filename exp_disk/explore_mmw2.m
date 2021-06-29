%% explor MMW disk model
sample=100;

mass=10.^(10);
md=0.1; %mstar/mv
fg=0;%[0.1 0.1];
beta=[0.5 0.2];
BT=0.1;
lambda=0.035;
xi=10.^(-1:4/99:3);


% check dependence on lambda
%lambda=lambda_prime(ones(sample,1));
%lambda=sort(lambda);
%sigma=zeros(length(mass),length(md),length(fg),length(beta),length(BT),length(xi),sample);
%rd=sigma;

%% build data base
tic
res=rscale_mmw_array(mass.*ones(sample,1),...
    'lambda',lambda,'md',md,'fg',fg(1),'beta',beta(1),...
    'BT',BT,'xi',xi,'noshow');
toc
rd=res.rd;
sigma=mass./(2*pi.*res.rd.^2);

tic
% res=rscale_mmw_array(mass.*ones(sample,1),...
%     'lambda',lambda,'md',md,'fg',fg(2),'beta',beta(2),...
%     'BT',BT,'xi',xi,'noshow');
% toc
% rd6=res.rd;
% sigma6=mass./(2*pi.*res.rd.^2);
% 


%     %% stellar + gas disk only
%
%     md=0.05:0.05:0.2; %mstar/mv
%     sigma=zeros(length(md),length(mass),sample);
%     rd=zeros(length(md),length(mass),sample);
%     tic
%     for i=1:length(md)
%         for j=1:length(mass)
%
%             res=rscale_mmw_array(mass(j).*ones(sample,1),'md',md(i),'lambda',lambda,'noshow');
%
%             rd(i,j,:)=res.rd;
%             sigma(i,j,:)=mass(j)./(2*pi.*res.rd.^2);
%
%
%             %sigma(i,j,1)=mean(sig);
%             %sigma(i,j,2)=std(sig);
%         end
%     end
%     toc
%
%
%
%
%
%
%
% end
% figure
% h=[];
% cc={'b' 'r' 'g' 'm' 'c'};
% lt={'.','s','o','+'};
% hold on
% for i=1:length(md)
%     for j=1:length(mass)
%         ll=sprintf('%s%s',cc{j},lt{i});
%         tag=sprintf('$M_s=%1.0e,M_v=%1.0e$',mass(j),mass(j)/md(i));
%         h(end+1)=loglog(lambda,squeeze(rd(i,j,:)),ll,'DisplayName',tag);
%
%     end
% end
% xlabelmine('$\lambda$')
% ylabelmine('$R_d\,[\mathrm{kpc}]$')
%
% hl=legend(h);
% set(hl,'Interpreter','latex','Location','NorthEastOutside','Fontsize',12)
%
%
%
%
%
% % cc={'b' 'r' 'g' 'm' 'c'};
% % hf=figure;
% % h=zeros(1,length(md));
% % hold on
% % for i=1:length(md)
% % h(i)=errorbar(mass,sigma(i,:,1),0.5.*sigma(i,:,2),sprintf('o%s',cc{i}),'DisplayName',sprintf('$m_d=%4.3f$',md(i)));
% % end
% %
% % grid
% % hl=legend(h);
% % set(hl,'Fontsize',12,'Interpreter','latex')
% % xlabelmine('$M_{star}\,[\mathrm{M_\odot}]$')
% % ylabelmine('$\Sigma_s\,[\mathrm{M_\odot/kpc^2}]$')