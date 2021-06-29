%% Sat
units;

% mSat=1e11;

%fc=0.1;
%fg=0.1;

% [rvSat,~,~,~]=calculate_virials('mv',mSat);
%
% pv=fc.*G/(8*pi).*(mSat.*Ms).^2./(rvSat.*Mpc).^4;
%% cluster
alfa=1; %0.5; %[0.5 1]; %[0.1 0.3 0.5 0.75 1.0];
%mSat=[1e10 1e11 1e12];

%[rvSat,~,~,~]=calculate_virials('mv',mSat);

%% get the shockedge values for the clusters
r0=0.1:0.1:3; % radius array in units of Rvir
aa='a1';
switch aa
    case 'a1'
        shEd=shockedge_a1;
    case 'a06'
        shEd=shockedge_a06;
    otherwise
        error('wronger pants')
end


eta1=[];
eta2=[];
h=[];
eta=0.01:0.01:50;
cc=[1 0 0; 0 0.7 0; 0 0 1]; % color definition array 
list=[101 102 103 104 105 106 107 3 5 6 7 9 10 11 14 24];
cnt=1;

%% run over selected clusters
for i= [1 ]; %  1:length(list)
    
    %cl=6;
    new_env(list(i),aa);
    
    %define the shock edge  
    sh1=shEd{i,4}./shEd{i,2};
    sh2=shEd{i,5}./shEd{i,2};
    sh=cat(2,sh1,sh2);
    shk=min(sh);
    
    % get virials
    mv(i)=get_mvir;
    rvC=get_rvir;
    vv=get_vvir;
    
    
    [mg,~,~]=read_Mass_Profiles(rvC);
    fc=mg./get_mvir; %define gas fraction of cluster
    
    
    rr=r0.*rvC; %radius in physical units 
    pad=7;
    [rog,~] = read_RHO_Profiles(rr); % density profile (replace with isothermal 
    %for j=1:length(alfa)
    
    pRPS=alfa.*(rog.*(Ms/Mpc.^3)).*(vv.*km).^2; % ram pressure of cluster 
    
    %% need to load catalog. 
     for j=1:1  %% go over different masses
%         switch j
%            case 1
%                cata=cata9_e2;
%            case 2
%                cata=cata10_e2;
%            case 3
%                cata=cata11_e2;
%         end
        %cata=cataSimclust100(j);
        
        %prepare stripping array 
        mstrip=zeros(length(cata.Ms),length(pRPS));
        mstrip2=zeros(length(cata.Ms),length(pRPS));
        
        
        hb = waitbar(0,'Processing catalog...');

             
        for k=1:length(cata.Ms)  % run over all galaxies 
            if mod(k,100)==0
            waitbar(k./length(cata.Ms),hb)
            end
            % find grav. force of galactic disk as a funciton of eta - the
            % position within the satellite 
            fd=disk_force_reduced(eta,'beta',cata.beta(k),'fg',cata.fg(k),...
                'BT',cata.BT(k));
            
            % create mhalo vs. eta from raw data
            mh0=interp1(cata.raw.rr(k,:)./cata.rd(k),cata.raw.mDM(i,:)./cata.Ms(k),eta,'pchip');
            fh=2.*mh0./eta.^2.*cata.beta(k).^2.*exp(-cata.beta(k).*eta);
            
            %total grav. force in the satellite 
            f1=(fd+fh).*pi.*G.*(cata.sigma(k).*(Ms/kpc^2)).^2.*cata.fg(k) ;
            mf1=exp_disk_mass(eta,cata.beta(k)); % mass fraction profile
            [f1Max, id]=max(f1); % max of force profile  - if the ram pressure is higher than this then you get total stripping
            
            
            for l=1:length(rr)
                                
                if pRPS(l)<=f1Max  % ram pressure less than max. grav, force 
                    ind=find(f1>pRPS(l),1,'last'); % indx of last stripped 
                    
                    % create small array for interpolation around stripping
                    % radius
                    l2=min(ind+pad,length(f1));
                    l1=max(1,ind-pad);
                    
                    ll=l1:l2;
                    
                    % define stripped mass  
                    mstrip(k,l)=interp1(f1(ll),mf1(ll),pRPS(l),'pchip');
                    
                else
                    mstrip(k,l)=0;  % total stripping 
                end
                
                f2=(fd+fh);
                [f2Max, id2]=max(f2);
                pRPS2=rps_factor_expdisk('alpha',alfa,'fc',fc,...
                    'fd',cata.fg(k),'rp',r0,'sigma',cata.sigma(k),...
                    'Mc',get_mvir);
                
                
                if pRPS2(l)<=f2Max
                    ind=find(f2>pRPS2(l),1,'last');
                    
                    l2=min(ind+pad,length(f2));
                    l1=max(1,ind-pad);
                    
                    ll=l1:l2;
                    
                    mstrip2(k,l)=interp1(f2(ll),mf1(ll),pRPS2(l),'pchip');
                    
                else
                    mstrip2(k,l)=0;
                end
                
            end
            
        end
        close(hb) 
        %mstrip=(1-mstrip).*100;
        mstrip=mstrip.*100;
        mstrip2=mstrip2.*100;
        
        
        cl(cnt).msat(j).mstrip=mstrip;
        cl(cnt).msat(j).mstrip2=mstrip2;
        cl(cnt).msat(j).msMean=mean(mstrip,1);
        cl(cnt).msat(j).msMax=max(mstrip,[],1);
        cl(cnt).msat(j).msMin=min(mstrip,[],1);
        cl(cnt).msat(j).msStd=std(mstrip,0,1);
        cl(cnt).msat(j).msMean2=mean(mstrip2,1);
        cl(cnt).msat(j).msMax2=max(mstrip2,[],1);
        cl(cnt).msat(j).msMin2=min(mstrip2,[],1);
        cl(cnt).msat(j).msStd2=std(mstrip2,0,1);
        cl(cnt).msat(j).logMass=log10(cata.Ms(1));
    end
    cl(cnt).shk=shk;
    cl(cnt).cl=list(i);
    cnt=cnt+1;
end

%% setup the ssfr proxy 
ssfr=zeros(size(mstrip));
sfr=zeros(size(mstrip));
for i=1:length(r0)
    mgas=cata.Ms.*cata.fg.*mstrip(:,i);
    ssfr(:,i)=ssfr_proxy(cata.Ms.*(1+cata.fb),mgas,3.*cata.rd);
    sfr(:,i)=ssfr(:,i).*cata.Ms.*(1+cata.fb);
end
ssfr=ssfr.*1e6; % convert to 1/Myr
sfr=sfr.*1e6;
units
fac=(Ms/Mpc^3)/(muMass*mp);
rog=rog.*fac;
%figure
c=[0 0 1;1 0 0; 0 0.7 0];
ms=[9 10 11];
%h=[];
for i=1:1
%     figure
%     h=[];
%     for j=1:1
%         l1=cl(i).msat(j).msMean+0.5.*cl(i).msat(j).msStd;
%         l2=cl(i).msat(j).msMean-0.5.*cl(i).msat(j).msStd;
%         l2=fliplr(l2);
%         ll=cat(2,l1,l2);
%         rl=cat(2,r0,fliplr(r0));
%         patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
%         hold on
%         h(j)=plot(r0,cl(i).msat(j).msMean,'color',c(j,:),'Linewidth',2,...
%             'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
%         %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
%         plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
%     end
%     
%     
%     hl=legend(h);
%     set(gca,'Fontsize',14,'box','on')
%     set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
%     
%     ylim([0 101])
%     % set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
%     grid
%     ylabelmine('Remaining Mass $[\%]$')
%     xlabelmine('$r/R_c$')
%     titlemine(sprintf('CL%s',num2str(cl(i).cl)))
%     
    
    figure 
    for j=1:1
        lines=prctile(mstrip,[25 50 70],1);
        l1=lines(1,:);
        l2=lines(3,:);
        l2=fliplr(l2);
        ll=cat(2,l1,l2);
        rl=cat(2,log10(rog),fliplr(log10(rog)));
        patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
        hold on
        h(j)=plot(log10(rog),lines(2,:),'color',c(j,:),'Linewidth',2,...
            'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
        %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
        %plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
    end
    
    
    %hl=legend(h);
    %set(gca,'Fontsize',14,'box','on')
    %set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
    
    %ylim([0 101])
    % set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
    grid
    ylabelmine('Remaining Mass $[\%]$')
    xlabelmine('$\log n_{\mathrm{gas}}\,[\mathrm{cm^{-3}}]$')%   r/R_c$')
    %titlemine(sprintf('CL%s',num2str(cl(i).cl)))
    titlemine('Mass Stripping vs. Density')
    
    
    figure 
    for j=1:1
        lines=prctile(mstrip,[25 50 70],1);
        l1=lines(1,:);
        l2=lines(3,:);
        l2=fliplr(l2);
        ll=cat(2,l1,l2);
        rl=cat(2,r0,fliplr(r0));
        patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
        hold on
        h(j)=plot(r0,lines(2,:),'color',c(j,:),'Linewidth',2,...
            'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
        %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
        %plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
    end
    
    grid
    ylabelmine('Remaining Mass $[\%]$')
    xlabelmine('$r/R_c$')
    titlemine('Mass Stripping vs. Radius')
    
     figure 
    for j=1:1
        lines=prctile(ssfr,[25 50 70],1);
        l1=lines(1,:);
        l2=lines(3,:);
        l2=fliplr(l2);
        ll=cat(2,l1,l2);
        rl=cat(2,log10(rog),fliplr(log10(rog)));
        patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
        hold on
        h(j)=plot(log10(rog),lines(2,:),'color',c(j,:),'Linewidth',2,...
            'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
        %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
        %plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
    end
    
    
    %hl=legend(h);
    %set(gca,'Fontsize',14,'box','on')
    %set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
    
    %ylim([0 101])
    % set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
    grid
    ylabelmine('SSFR $[\mathrm{Myr^{-1}}]$')
    xlabelmine('$\log n_{\mathrm{gas}}\,[\mathrm{cm^{-3}}]$')
    titlemine('SSFR vs. Density')
    
    figure 
    for j=1:1
        lines=prctile(ssfr,[25 50 70],1);
        l1=lines(1,:);
        l2=lines(3,:);
        l2=fliplr(l2);
        ll=cat(2,l1,l2);
        rl=cat(2,r0,fliplr(r0));
        patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
        hold on
        h(j)=plot(r0,lines(2,:),'color',c(j,:),'Linewidth',2,...
            'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
        %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
        %plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
    end
    
    
    %hl=legend(h);
    %set(gca,'Fontsize',14,'box','on')
    %set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
    
    %ylim([0 101])
    % set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
    grid
     ylabelmine('SSFR $[\mathrm{Myr^{-1}}]$')
    xlabelmine('$r/R_c$')
    titlemine('SSFR vs. Radius')
    
    figure 
    for j=1:1
        lines=prctile(sfr,[25 50 70],1);
        l1=lines(1,:);
        l2=lines(3,:);
        l2=fliplr(l2);
        ll=cat(2,l1,l2);
        rl=cat(2,log10(rog),fliplr(log10(rog)));
        patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
        hold on
        h(j)=plot(log10(rog),lines(2,:),'color',c(j,:),'Linewidth',2,...
            'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
        %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
        %plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
    end
    
    
    %hl=legend(h);
    %set(gca,'Fontsize',14,'box','on')
    %set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
    
    %ylim([0 101])
    % set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
    grid
    ylabelmine('SFR $[\mathrm{M_{\odot}\,Myr^{-1}}]$')
    xlabelmine('$\log n_{\mathrm{gas}}\,[\mathrm{cm^{-3}}]$')
    titlemine('SFR vs. Density')
    
    figure 
    for j=1:1
        lines=prctile(sfr,[25 50 70],1);
        l1=lines(1,:);
        l2=lines(3,:);
        l2=fliplr(l2);
        ll=cat(2,l1,l2);
        rl=cat(2,r0,fliplr(r0));
        patch(rl,ll,c(j,:),'faceAlpha',0.2,'EdgeColor','none')
        hold on
        h(j)=plot(r0,lines(2,:),'color',c(j,:),'Linewidth',2,...
            'DisplayName',sprintf('$M_s=10^{%s}\\,\\mathrm{M_\\odot}$',num2str(ms(j))));
        %'DisplayName',sprintf('CL%s',num2str(cl(i).cl)));
        %plot(cl(i).shk.*[1 1],[0 100],'--k','linewidth',3)
    end
    
    
    %hl=legend(h);
    %set(gca,'Fontsize',14,'box','on')
    %set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast','box','on')
    
    %ylim([0 101])
    % set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
    grid
     ylabelmine('SFR $[\mathrm{M_{\odot}\,Myr^{-1}}]$')
    xlabelmine('$r/R_c$')
    titlemine('SFR vs. Radius')
    
    
    
end
%eta1(i,k,j)=interp1(r0,etaS,1);
%eta2(i,k,j)=interp1(r0,etaS,2);

%             switch i
%                 case 1
%                     ci=1;
%                 case 8
%                     ci=2;
%                 case 16
%                     ci=3;
%                 otherwise
%                    error('wrong pants')
%             end
%
%
%
%             h(end+1)=plot(r0,etaS,'linewidth',2,'color',cc(ci,:),...
%                 'DisplayName',sprintf('CL%s,%s',num2str(list(i)),mk_mvir_string(mv(i))));
%             hold on
%             h(end+1)=plot(r0,etaIso,'--','linewidth',1.5,'color',cc(ci,:),...
%                 'DisplayName',sprintf('CL%s, Isothermal',num2str(list(i))));
%

% hl=legend(h([1 3 5]));
% set(hl,'Fontsize',14','Interpreter','latex','Location','NorthWest')
% ylim([0 0.4])
% set(gca,'Fontsize',14,'box','on','Ytick',[0:0.05:0.4])
% grid
% xlabelmine('$r_p/R_c$')
% ylabelmine('$r_s/R_s$')
%
%
% figure
% h=[];
% h(1)=semilogx(mv,1-eta1(:,2,1),'sk','markerFaceColor','b','markerSize',12,...
%     'DisplayName','$1 R_{\mathrm{vir}}$');
% hold on
% h(2)=semilogx(mv,1-eta2(:,2,1),'ok','markerFaceColor','r','markerSize',12,...
%     'DisplayName','$2R_{\mathrm{vir}}$');
% xlim([5e13 1e15])
% grid
%
% hl=legend(h);
%  set(hl,'Fontsize',14','Interpreter','latex','Location','SouthEast')
% % ylim([0 0.4])
%  set(gca,'Fontsize',14,'box','on')
%
%  xlabelmine('$M_{\mathrm{vir}}\,[\mathrm{M_\odot}]$')
%  ylabelmine('Mass Stripped [\%]')
%
%
%
%
