%% Sat
units;

% mSat=1e11;

%fc=0.1;
fg=0.1;

% [rvSat,~,~,~]=calculate_virials('mv',mSat);
%
% pv=fc.*G/(8*pi).*(mSat.*Ms).^2./(rvSat.*Mpc).^4;
%% cluster
alfa=[0.5 1]; %[0.1 0.3 0.5 0.75 1.0];
mSat=[1e10 1e11 1e12];

[rvSat,~,~,~]=calculate_virials('mv',mSat);
r0=0.5:0.01:3;

aa='a06';
switch aa
    case 'a1'
        shEd=shockedge_a1;
        shkM1=shockedgeMask1_a1;
        shkM2=shockedgeMask2_a1;
        b=0;
    case 'a06'
        shEd=shockedge_a06;
        shkM1=shockedgeMask1_a06;
        shkM2=shockedgeMask2_a06;
        b=1;
    otherwise
        error('wronger pants')
end


hf=figure;
eta1=[];
eta2=[];
h=[];

cc=[1 0 0; 0 0.7 0; 0 0 1];
cc=brewermap(3,'Set1');
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
for i= 1:length(list) %[1 8 15];%  
    
    sh1=shEd{i,4}./shEd{i,2};
    sh2=shEd{i,5}./shEd{i,2};
    
    
    sh=cat(2,sh1(shkM1(i,:)),sh2(shkM2(i,:)));
    shkMn=min(sh);
    shkMx=max(sh);
    
    %cl=6;
    new_env(list(i),aa);
    global zred
    
    
    mv(i)=get_mvir;
    rvC=get_rvir;
    vv=get_vvir;
    
    [mg,~,~]=read_Mass_Profiles(shkMn.*rvC);
    fc=fg; %mg./get_mvir;
      
    
    rr=r0.*rvC;
    
    [rog,~] = read_RHO_Profiles(rr);
%      figure(10+i+b)
%      loglog(r0,rog)
   
    
    ls=0.01:0.0001:2; 
    for j=   1:1 %length(alfa)
        pRPS=alfa(j).*(rog.*(Ms/Mpc.^3)).*(vv.*km).^2;
        
        cvC=cvir_Mvir(mv(i),zred);
        pRPSN0=alfa(j).*fg.*G/(4*pi).*(mv(i).*Ms).^2./(rvC.*Mpc).^4;
        pRPSN=pRPSN0./(nfwA(1,cvC).*r0.*(cvC^-1+r0).^2);
        for k=2:2 %length(mSat)
             
            cvSat=cvir_Mvir(mSat(k),zred);
            as=nfwA(1,cvSat);
           
            %pv=fg.*G/(4*pi).*(mSat(k).*Ms).^2./(rvSat(k).*Mpc).^4;
            
            pv1=fg.*G/(4*pi).*(mSat(k).*Ms).^2./(rvSat(k).*Mpc).^4;
            pv2=nfwA(ls,cvSat)./(as.*ls.*(cvSat.^-1+ls)).^2;
            
            pv=pv1.*pv2;
            
            
            
            
            
            for m=1:length(r0)
                lStrip(m)=interp1(pv,ls,pRPS(m),'PCHIP');
                
                lStrip2(m)=interp1(pv,ls,pRPSN(m),'PCHIP');
                
            end
            mStrip=nfwA(lStrip,cvSat)./nfwA(1,cvSat);            
            mStrip2=nfwA(lStrip2,cvSat)./nfwA(1,cvSat);            
            
            
            
            
            %etaS=sqrt(pv./pRPS);
            
            %etaIso=1/sqrt(alfa(j)).*sqrt(fg./fc).*(mSat(k)./mv(i)).^(1/3).*r0;
            
            
            eta1(i,k,j)=interp1(r0,mStrip,1);
            eta2(i,k,j)=interp1(r0,mStrip,2);
            
            switch i
                case 1 
                    ci=1;
                case 8
                    ci=2;
                case 15
                    ci=3;
                otherwise
                   %error('wrong pants')
            end
                   
%             
%             figure(hf)        
%             h(end+1)=plot(r0,mStrip,'linewidth',2,'color',cc(ci,:),...
%                 'DisplayName',sprintf('CL%s,%s',num2str(list(i)),mk_mvir_string(mv(i))));
%             hold on
%             plot(r0,mStrip2,'--','linewidth',1.5,'color',cc(ci,:),...
%                 'DisplayName',sprintf('CL%s, Isothermal',num2str(list(i))));
%             plot(shkMn.*[1 1],[0 1],':','linewidth',3,'color',cc(ci,:));
%             %plot(shkMx.*[1 1],[0 1],'.-','linewidth',3,'color',cc(ci,:));
%             
        end
    end
end
%hl=legend(h([1 4 7]));
% hl=legend(h);
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% xlim([0.5 3])
% ylim([0 0.9])
% set(gca,'Fontsize',14,'box','on','Ytick',0:0.1:0.9)
% grid 
% xlabelmine('$r_p/R_c$')
% ylabelmine('$\ell_s/R_{sat}$')
% ylabelmine('$M(l<l_{\mathrm{s}})/M_{\mathrm{sat}}$')


figure
h=[];
h(1)=semilogx(mv,100.*(1-eta1(:,2,1)),'sk','markerFaceColor','b','markerSize',12,...
    'linewidth',1.5,'DisplayName','$1 R_{\mathrm{vir}}$');
hold on
h(2)=semilogx(mv,100.*(1-eta2(:,2,1)),'ok','markerFaceColor','r','markerSize',12,...
    'linewidth',1.5,'DisplayName','$2R_{\mathrm{vir}}$');
switch aa
    case 'a1'
    xlim([5e13 3e15])
    case 'a06'
        xlim([5e13 1e15])
    otherwise
        error('wtf')
end
   
grid 

hl=legend(h);
 set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast')
% ylim([0 0.4])
 set(gca,'Fontsize',14,'box','on')

 xlabelmine('$M_{\mathrm{vir}}\,[\mathrm{M_\odot}]$')
 ylabelmine('Mass Stripped [\%]')




