%output_dir='C:/Users/owner/Ubuntu One/cluster/printout';

%% basic plot, alpha=1

% basic values
zred=0;
mc=1e15;
ms=[1e12 1e13 1e14] ;

cvfac=[0.9 1 1.1];
cc=cvfac.*cvir_Mvir(mc,zred);


ac=nfwA(1,cc);


fg=0.1;
fc=0.1;
alfa=0.5;

rp=0.05:0.001:3.1; %in units of cluser rv


ls=0.01:0.0001:2;
% get lhs of stripping equation for each satellite
%lhs=zeros(size(ls));

lStrip=zeros(size(rp));
msN=zeros(length(rp),length(ms),length(cvfac));



for j=1:length(ms);
    %cs=fliplr(cvfac).*cvir_Mvir(ms(j),zred);
    cs=cvfac.*cvir_Mvir(ms(j),zred);
    as=nfwA(1,cs);
    
    for k=1:length(cvfac)
        
        
        lhs=(ls.^2.*(cs(k).^-1+ls).^2)./nfwA(ls,cs(k)); % radius in satellite
        
        
        for i=1:length(rp)
            
            rhs=(fg*(1+0)/fc/alfa).*(ms(j)./mc).^(2/3).*(ac(k)./as(k).^2).*rp(i).*(cc(k)^-1+rp(i)).^2;%./nfwA(rc,cc);
            
            ff=lhs-rhs;
            
            lStrip(i)=interp1(ff,ls,0,'PCHIP');
            %lStrip(i,j)=interp1(ff,ls,0,'PCHIP');
            %ii=find(ff>0,1,'first');
            %lStrip(i,j)=ls(ii);
            
        end
        
        msN(:,j,k)=nfwA(lStrip,cs(k))./nfwA(1,cs(k));
    end
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
    %hN(j)=plot(rp,msN(:,j,1),':','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    %hN(j)=plot(rp,msN(:,j,3),':','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
    hI(j)=plot(rp,msI(:,j),'--','color',cc(j,:),'linewidth',2);%,'DisplayName',dnTag);
    hold on
end

for j=1:length(ms)
    dnTag=sprintf('$M_{\\mathrm{C}},m_{\\mathrm{sat}}=10^{15},10^{%s}$',num2str(floor(log10(ms(j)))));
    
   
    hold on
    %plot(rp,msN(:,j,1),'-','color',cc(j,:),'linewidth',1.5,'DisplayName',dnTag);
    %plot(rp,msN(:,j,3),'-','color',cc(j,:),'linewidth',1.5,'DisplayName',dnTag);
    xx=cat(2,rp,fliplr(rp));
    yy=cat(2,msN(:,j,1)',fliplr(msN(:,j,3)'));
    patch(xx,yy,cc(j,:),'facealpha',0.3,'linestyle','none','DisplayName',dnTag);
    hold on
     hN(j)=plot(rp,msN(:,j,2),'-','color',cc(j,:),'linewidth',2,'DisplayName',dnTag);
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


