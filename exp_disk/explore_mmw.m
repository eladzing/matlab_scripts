%% explor MMW disk model
if buildFlag
    sample=30;
    
    mass=10.^(8:1:11);
    md=0.05:0.05:0.15; %mstar/mv
    fg=0.0:0.1:0.2;
    beta=[0.5 1.0 2.0];
    BT=0:0.1:0.2;
    xi=[1 5];
    
    
    % check dependence on lambda
    lambda=lambda_prime(ones(sample,1));
    
    sigma=zeros(length(mass),length(md),length(fg),length(beta),length(BT),length(xi),sample);
    rd=sigma;
    
    %% build data base
    
    tIter=0.75;
    tEst=tIter*numel(rd);
    
    hr=floor(tEst/3600);
    mins=floor((tEst-hr*3600)/60);
    sec=ceil(tEst-hr*3600-mins*60);
    
    hwb=waitbar(0,sprintf('remaning time: %i hours, %i minutes %i seconds',hr,mins,sec),...
        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(hwb,'canceling',0)
    
    tAll=tic;
    for i=1:length(mass)
        tStart=tic;
        for j=1:length(md)
            for k=1:length(fg)
                for l=1:length(beta)
                    for m=1:length(BT)
                        for n=1:length(xi)
                            % Check for Cancel button press
                            if getappdata(hwb,'canceling')
                                break
                            end
                            res=rscale_mmw_array(mass(i).*ones(sample,1),'lambda',lambda,...
                                'md',md(j),'fg',fg(k),'beta',beta(l),...
                                'BT',BT(m),'xi',xi(n),'noshow');
                            
                            rd(i,j,k,l,m,n,:)=res.rd;
                            sigma(i,j,k,l,m,n,:)=mass(j)./(2*pi.*res.rd.^2);
                            
                        end
                    end
                end
            end
        end
        tElapsed=toc(tStart);
        tEstimate=(length(mass)-i).*tElapsed;
        hr=floor(tEstimate/3600);
        mins=floor((tEstimate-hr*3600)/60);
        sec=ceil(tEstimate-hr*3600-mins*60);
        
        waitbar(i/length(mass),hwb,sprintf('remaning time: %i hours, %i minutes %i seconds',hr,mins,sec ));
        
    end
    toc(tAll)
    delete(hwb)
end
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