global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/wetzel.mat'])
RomC_data;

load('/home/zinger/workProjects/matlab_scripts/rpsProject/matFiles/qFracFullOrb.mat')

betaRange=[1/5 1/3 1/2 1/1.5 4/5 1.0 Inf]; 
fgsRange=[0.05 0.1 0.2 0.3 0.4 0.5 Inf];
    



lw=1.5;
%betaTag='$\beta';
beta={'1/5' '1/3' '1/2' '2/3' '4/5'...
    '1' 'Rand'};
beta2={'0.2' '0.33' '0.5' '0.66' '0.8' '1'...
    'Rand'};
betaPrint={'02' '033' '05' '066' '08' '1'...
    'Rand'};
%fgsTag='$f_\mathrm{gs}';
fgs={'0.05' '0.1' '0.2' ...
    '0.3' '0.4' '0.5' ...
    'Rand'};
fgsPrint={'005' '01' '02' ...
    '03' '04' '05' ...
    'Rand'};

colo=brewermap(8,'Set1');
colo=colo([1:5 7],:);

%% keep fgs constant, explore beta 
for i=1:7
    figure

    h=[];
    for j=1:7
        if j~=7
            bTag=['$\beta=' beta{j} '$'];
            col=colo(j,:);
            st='-';
        else
            bTag=[beta{j} ' $\beta$'];
            st='--';
            col='k';
        end
        bCen=qFrac(i,j).bCen;
        qf=qFrac(i,j).qfProj./qFrac(i,j).cntProj;
        h(j)=plot(bCen,qf,st,'color',col,...
            'linewidth',lw,'Displayname',bTag);
        if j==1; hold on;  end
        % h(2)=plot(bCen,qf2./cnt2,'rx-',...
        %        'linewidth',lw,'Displayname','projected');
    end
    
    h(end+1)=plot(0.68.*Wetzel12_fig15(2,:),Wetzel12_fig15(1,:),'sk',...
        'DisplayName','Wetzel+12','markersize',8);
%     h(end+1)=errorbar(fig15.rdist,fig15.qf,fig15.errM,fig15.errP,'ok',...
%      'DisplayName','RomulusC,$>10^{9.7}$');
%     h(end+1)=errorbar(rdistBLow,fig15.qfLow,fig15.errMlow,fig15.errPlow,'ok',...
%      'DisplayName','RomulusC,$<10^{9.7}$');

    ylim([0 1])
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',12,'Location','southwest','NumColumns',3);
    
    
    if i~=7
        tTag=['$f_\mathrm{gs}=' fgs{i} '$'];
    else
        tTag=[fgs{i} ' $f_\mathrm{gs}$'];
    end
    
    
    
    xlabelmine('$r/R_{\mathrm{200,c}}$',16);
    ylabelmine('quenched fraction',16);
    xlim([0 1.05])
    %titlemine(sprintf('$f_{\\mathrm{gs}}=%s$',fgs{i}));
    titlemine(tTag);
    set(gca,'fontsize',14)
    fname=['qFrac_OrbEnsembleFull_surv_beta_fgs' fgsPrint{i}];
    printout_fig(gcf,fname,'v')
end

%% keep beta constant, explore fgs 
for i=1:7
    figure
    
    h=[];
    for j=1:7
        if j~=7
            bTag=['$f_\mathrm{gs}=' fgs{j} '$'];
            col=colo(j,:);
            st='-';
        else
            bTag=[fgs{j} ' $f_\mathrm{gs}$'];
            st='--';
            col='k';
        end
        
        
        bCen=qFrac(j,i).bCen;
        qf=qFrac(j,i).qfProj./qFrac(j,i).cntProj;
        h(j)=plot(bCen,qf,st,'color',col,...
            'linewidth',lw,'Displayname',bTag);
        if j==1; hold on;  end
        % h(2)=plot(bCen,qf2./cnt2,'rx-',...
        %        'linewidth',lw,'Displayname','projected');
    end
    
    h(end+1)=plot(0.68.*Wetzel12_fig15(2,:),Wetzel12_fig15(1,:),'sk',...
        'DisplayName','Wetzel+12','markersize',8);
%       h(end+1)=errorbar(fig15.rdist,fig15.qf,fig15.errM,fig15.errP,'ok',...
%      'DisplayName','RomulusC,$>10^{9.7}$');
%     h(end+1)=errorbar(rdistBLow,fig15.qfLow,fig15.errMlow,fig15.errPlow,'ok',...
%      'DisplayName','RomulusC,$<10^{9.7}$');
    
    xlim([0 1.05])
    ylim([0 1])
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',12,'Location','southwest','NumColumns',3);
    
    if i~=7
        tTag=['$\beta=' beta{i} '$'];
    else
        tTag=[beta{i} ' $\beta$'];
    end
    
    
    xlabelmine('$r/R_{\mathrm{200,c}}$',16);
    ylabelmine('quenched fraction',16);
    titlemine(tTag);
    set(gca,'fontsize',14)
    fname=['qFrac_OrbEnsembleFull_surv_fgs_beta' betaPrint{i}];
    printout_fig(gcf,fname,'v')
end



