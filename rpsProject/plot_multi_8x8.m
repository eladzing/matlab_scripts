global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/wetzel.mat'])


lw=1.5;
%betaTag='$\beta';
beta={'1/2.5' '1/2' '2/3' '1' '1.5'...
    '2' '2.5' 'Rand'};
beta2={'0.4' '0.5' '0.66' '1' '1.5'...
    '2' '2.5' 'Rand'};
%fgsTag='$f_\mathrm{gs}';
fgs={'0.05' '0.2' '0.35' ...
    '0.5' '0.65' '0.8' ...
    '0.95' 'Rand'};

colo=brewermap(8,'Set1');

%% keep fgs constant, explore beta 
for i=1:8
    figure

    h=[];
    for j=1:8
        if j~=8
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
    
    h(end+1)=plot(Wetzel12_fig15(2,:),Wetzel12_fig15(1,:),'sk',...
        'DisplayName','Wetzel+12','markersize',8);
    
    ylim([0 1])
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',12,'Location','southwest','NumColumns',3);
    
    
    if i~=8
        tTag=['$f_\mathrm{gs}=' fgs{i} '$'];
    else
        tTag=[fgs{i} ' $f_\mathrm{gs}$'];
    end
    
    
    
    xlabelmine('$r/R_{\mathrm{200,m}}$',16);
    ylabelmine('quenched fraction',16);
    %titlemine(sprintf('$f_{\\mathrm{gs}}=%s$',fgs{i}));
    titlemine(tTag);
    set(gca,'fontsize',14)
    fname=['qFrac_radOrb_surv_beta_fgs' fgs{i}];
    printout_fig(gcf,fname,'v')
end

%% keep beta constant, explore fgs 
for i=1:8
    figure
    
    h=[];
    for j=1:8
        if j~=8
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
    
    h(end+1)=plot(Wetzel12_fig15(2,:),Wetzel12_fig15(1,:),'sk',...
        'DisplayName','Wetzel+12','markersize',8);
    
    ylim([0 1])
    grid
    
    hl=legend(h);
    set(hl,'Interpreter','latex','fontsize',12,'Location','southwest','NumColumns',3);
    
    if i~=8
        tTag=['$\beta=' beta{i} '$'];
    else
        tTag=[beta{i} ' $\beta$'];
    end
    
    
    xlabelmine('$r/R_{\mathrm{200,m}}$',16);
    ylabelmine('quenched fraction',16);
    titlemine(tTag);
    set(gca,'fontsize',14)
    fname=['qFrac_radOrb_surv_fgs_beta' beta2{i}];
    printout_fig(gcf,fname,'v')
end



