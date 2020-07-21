 %output_dir='C:/Users/owner/Ubuntu One/cluster/printout';

%% basic plot, alpha=1   

% basic values 
zred=0;
mc=1e15;
ms=[1e10 1e11 1e12];

cc=1.1.*cvir_Mvir(mc,zred);
cs=0.9.*cvir_Mvir(ms,zred);

ac=nfwA(1,cc);
as=nfwA(1,cs);

fg=0.1;
fc=0.1;
alfa=0.5;

rp=0.1:0.001:3.1; %in units of cluser rv


ls=0.01:0.0001:2; 
% get lhs of stripping equation for each satellite
%lhs=zeros(size(ls));
lStrip=zeros(length(rp),length(ms));
msN=zeros(length(rp),length(ms));

for j=1:length(ms);
    lhs=(ls.^2.*(cs(j).^-1+ls).^2)./nfwA(ls,cs(j)); % radius in satellite
    
    
    for i=1:length(rp)
                
        rhs=(fg*(1+0)/fc/alfa).*(ms(j)./mc).^(2/3).*(ac./as(j).^2).*rp(i).*(cc^-1+rp(i)).^2;%./nfwA(rc,cc);
        
        ff=lhs-rhs;
        
        lStrip(i,j)=interp1(ff,ls,0,'PCHIP');
        %ii=find(ff>0,1,'first');
        %lStrip(i,j)=ls(ii);
        
    end
    
    msN(:,j)=nfwA(lStrip(:,j),cs(j))./nfwA(1,cs(j));
end


%% include Isothermal Model 

massRatio=ms./mc; %  [1e-5 1e-4 1e-3]; % virial mass ratio
msI=zeros(length(rp),length(massRatio));
for i=1:length(massRatio)
      msI(:,i)=massRatio(i).^(1/3).*rp./sqrt(alfa);
end


%% plot 
cc=brewermap(8,'Set1');
figure
hN=[];
hI=[];
for j=1:length(ms)
    dnTag=sprintf('$m_{\\mathrm{sat}}=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(floor(log10(ms(j)))));
    hN(j)=plot(rp,msN(:,j),'-','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    hold on
    hI(j)=plot(rp,msI(:,j),'--','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    
end
    
grid
xlim([0.5 3])
ylim([1e-2 0.8])

%hl=legend(cat(2,hN,hI));
hl=legend(hN);

set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);

set(gca,'Fontsize',12);
xlabelmine('$r_p/R_c$')
ylabelmine('$M(l<l_{\mathrm{s}})/M_{\mathrm{sat}}$')

    %plot(rp,rs,'Linewidth',2)

% grid
% ylim([1e-2 0.5])
% hl=legend('$m_{sat}=10^{9}\,\mathrm{M_\odot}$','$m_{sat}=10^{10}\,\mathrm{M_\odot}$','$m_{sat}=10^{11}\,\mathrm{M_\odot}$','$m_{sat}=10^{12}\,\mathrm{M_\odot}$');%,'$m_{sat}=10^{13}\,\mathrm{M_\odot}$');
% set(hl,'Interpreter','latex','Location','NorthWest','Fontsize',14);
% set(gca,'Fontsize',12);
% xlabelmine('$r_p/R_c$')
% ylabelmine('$r_s/R_s$')
% titlemine('Stripping radius in a $10^{14}\,\mathrm{M_\odot}$ Cluster ($\alpha=1$)')
% name='%s/rps_toy_m14.%s';
% %exportfig(gcf,sprintf(name,output_dir,'png'),'format','png');
% %exportfig(gcf,sprintf(name,output_dir,'eps'));


