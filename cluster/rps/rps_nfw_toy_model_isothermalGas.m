 

%% plot for exploring mu and alpha dependance 
% basic values  

zred=0;
mc=1e15;

mu=-5:0.01:-3; % virial mass ratio
ms=10.^mu.*mc;


cc=cvir_Mvir(mc,zred);
cs=cvir_Mvir(ms,zred);

ac=nfwA(1,cc);
as=nfwA(1,cs);

fg=0.1;
fc=0.1;
alfa=0.5;

xx=0.0001:0.0001:1;


alpha=0.01:0.01:1.0;

rp=2.0;
ls=0.01:0.0001:2; 
mstrip=zeros(length(alpha),length(ms));
for j=1:length(ms)
    
    
    bf=(trapz(xx,xx.^2.*exp(2./as(j).*log(cs(j).*xx+1)./xx))).^-1;
    
    msat=bf.*cumtrapz( xx, xx.^2.*exp(2./as(j).*log(cs(j).*xx+1)./xx));
    
    rhs0=(ms(j)./mc).^(2/3).*(ac./as(j)).*bf.*exp(2/as(j).*log(1+cs(j).*ls)./ls); %%(ls.^2.*(cs(k).^-1+ls).^2)./nfwA(ls,cs(k)); % radius in satellite
    lhs=1./( rp.*(cc^-1+rp).^2);%./nfwA(rc,cc);
    
    
    %%lhs=(ls.^2.*(cs(j).^-1+ls).^2)./nfwA(ls,cs(j)); % radius in satellite
    
    
    for i=1:length(alpha)
        
        rhs=rhs0.*(fg*(1+0)/fc/alfa);
        
        ff=rhs./lhs;
        
        lStrip=interp1(ff,ls,1,'PCHIP');
        %rhs=(fg*(1+0)/fc/alpha(i)).*(ms(j)./mc).^(2/3).*(ac./as(j).^2).*rp.*(cc^-1+rp).^2;%./nfwA(rc,cc);
        %  ff=lhs-rhs;
        
        %%ii=find(ff>0,1,'first');
        %%lStrip=interp1(ff,ls,0,'PCHIP');%      ls(ii);
        
        %if ~isempty(lStrip)
            
        mstrip(i,j)=interp1(xx,msat,lStrip,'PCHIP'); 
        %mstrip(i,j)=nfwA(lStrip,cs(j))./nfwA(1,cs(j));
        %else
        %    disp('*')
        %    mstrip(i,j)=1;
        %end
        
        
    end
end
figure
map = brewermap(256,'*Spectral');
colormap(map)
imagesc(log10(ms),alpha,(mstrip),[0 1]);
set(gca,'Ydir','normal','FontSize',14,'Ytick',0:0.1:1);%,'Xtick',-5:0.5:-2);
bar=colorbar;

barTitle(bar,'$M(l<l_{\mathrm{s}})/M_{\mathrm{sat}}$');
%set(get(bar,'Title'),'String',,'Fontsize',12,'Interpreter','latex');
%set(bar,'FontSize',12,'Ytick',0:10:100)
set(bar,'FontSize',12,'Ytick',0:0.1:1)
xlabelmine('$\log(M_{\mathrm{sat}})\,[\mathrm{M_{\odot}}]$')
ylabelmine('$\varepsilon^{\prime}$')
%titlemine('Mass strippes at $r_p=1\,R_{\mathrm{vir}}$')
