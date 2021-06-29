%% explor MMW disk model
sample=1;

mass=10.^(8);
md=0.15; %:-0.05:0.05; %mstar/mv
lambda=0.02:-0.001:0.001;
fg=0.5;
beta=0.5;
BT=0;
xi=100;


% check dependence on lambda
%lambda=lambda_prime(ones(sample,1));
%lambda=sort(lambda);
%sigma=zeros(length(mass),length(md),length(fg),length(beta),length(BT),length(xi),sample);
%rd=sigma;

%% build data base
tElap=zeros(length(md),length(lambda));
for j=1:length(md)
    for i=1:length(lambda)
        fprintf('md=%5.3f, lambda=%5.3f  \n',md(j),lambda(i));
        tic
        res=rscale_mmw_array(mass.*ones(sample,1),...
            'lambda',lambda(i),'md',md(j),'fg',fg,'beta',beta,...
            'BT',BT,'xi',xi,'noshow');
        tElap(j,i)=toc;
    end
    
end

%rd=res.rd;
%sigma=mass./(2*pi.*res.rd.^2);
%toc


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