%% batch plot properties
printFlag=true;

units;
sims={'TNG50','TNG100'};


xfields={'galStellarMass','hostM200c','rpos','vrad','vel','galGasMass','massRatio','radiality'};
xlab={'log Stellar Mass $[\mathrm{M_\odot}]$',...
    '$\log M_\mathrm{200,c}\, [\mathrm{M_\odot}]$',...
    'Radial Position $[R_\mathrm{200,c}]$',...
    'Radial Velocity $[V_\mathrm{200,c}]$',... 
    'log Velocity $[V_\mathrm{200,c}]$',... 
    'log Gas Mass $[\mathrm{M_\odot}]$',...
    '$\log M_\mathrm{sat}/M_\mathrm{host}$',...
    '$v_\mathrm{rad}/|\vec{v}|$'};

xlog=false(size(xfields));
xlog([1 2  5 6 7 ])=true;

yfields={'hostM200c','rpos','vrad','vel','galSFR','galGasMass','gasMass',...
    'galBHMass','massRatio','radiality'};
ylab={'$\log M_\mathrm{200,c}\, [\mathrm{M_\odot}]$',...
    'Radial Position $[R_\mathrm{200,c}]$',...
    'Radial Velocity $[V_\mathrm{200,c}]$',...
    'log Velocity $[V_\mathrm{200,c}]$',...
    'log sSFR $[\mathrm{yr^{-1}}]$',...
    'log Gas Mass $[\mathrm{M_\odot}]$',...
    'log CGM Mass $[\mathrm{M_\odot}]$',...
    'log BH Mass $[\mathrm{M_\odot}]$',...
    '$\log M_\mathrm{sat}/M_\mathrm{host}$',...
    '$v_\mathrm{rad}/|\vec{v}|$'};

skip=false(length(xfields),length(yfields));
skip(:,1)=true;
skip(4,2:4)=true;
skip(5,2)=true;
skip(6,1:4)=true;
skip(7,1:4)=true;
skip(8,[1:4 6 9])=true;



ylog=false(size(yfields));
ylog([1 4 5 6 7 8 9])=true;
% some perliminaries
mv=galProps.hostM200c;
rv=galProps.hostR200c;
vvir=sqrt(Units.G.*(mv.*Units.Ms)./(rv.*Units.kpc))./Units.km; %in km/sec

sSFR=galProps.galSFR./galProps.galStellarMass+1e-13;
vv=(sqrt(sum(galProps.vel.^2,1))./vvir);
vr=galProps.vrad./vvir;
rp=galProps.rpos./rv;
cgmMass=galProps.gasMass-galProps.galGasMass;
massRatio=galProps.galStellarMass./mv;
radiality=galProps.vrad./vv;


%% plotting defaults 
cmap=brewermap(256,'Reds');
jfCol=[0.1,0.1,0.89];
njfCol=[0.89,0.1,0.11];
cols=cat(1,njfCol,jfCol);
otherCol=brewermap(8,'Set1');

outdir='jellyfish';


for k=1:length(sims)
    
    switch sims{k}
        case 'TNG50'
            simMask=mask50;
        case 'TNG100'
            simMask=mask100;
    end
    
    if k==1
        startwithI=3;
    else
        startwithI=1;
    end
    
    for i=startwithI:length(xfields)
        
        switch xfields{i}
            case 'rpos'
                xx0=rp;
            case 'vrad'
                xx0=vr;
            case 'vel'
                xx0=vv;
            case 'massRatio'
                xx0=massRatio;
            case 'radiality'
                xx0=radiality;
            otherwise
                xx0=galProps.(xfields{i});
        end
        if xlog(i)
            xx0=log10(xx0);
        end
            
        xx=xx0(simMask);
        xxJF=xx0(simMask & maskJF);
        xxNJF=xx0(simMask & ~maskJF);
        
        if i==1
            startwith=1;
        else
            startwith=1;
        end
        
        for j=startwith:length(yfields)
            if strcmp(xfields{i},yfields{j})
                continue
            end
            
            if skip(i,j)
                continue;
            end
            
            switch yfields{j}
                case 'rpos'
                    yy0=rp;
                case 'vrad'
                    yy0=vr;
                case 'vel'
                    yy0=vv;
                case 'massRatio'
                    yy0=massRatio;
                case 'galSFR'
                    yy0=sSFR;
                case 'gasMass'
                    yy0=cgmMass;
                    yy0(yy0<=0)=1e4;
                case 'radiality'
                    yy0=radiality;
                otherwise
                    yy0=galProps.(yfields{i});
            end
           
            if ylog(j)
                
                yy0=log10(yy0);
            end
            
            yy=yy0(simMask);
            yyJF=yy0(simMask & maskJF);
            yyNJF=yy0(simMask & ~maskJF);
            
                        
            %% prepare data 
            xl0=[min(xx) max(xx)].*[0.99 1.01];
            yl0=[min(yy) max(yy)].*[0.99 1.01];
            
            [bird, binsize, xl,yl]= histogram2d(xxNJF,yyNJF,ones(size(xxNJF)),...
                'xlim',xl0,'ylim',yl0,'len',50);
            
            filt=fspecial('disk',15);
            popCont=plot_population_contour(xxJF,yyJF,'smooth',filt,'noplot');
            
            jfscore=maskJF(simMask);
            
            %% identify points beyond contour
            
            
            lx=max(diff(popCont.xx));ly=max(diff(popCont.yy));
            x0=popCont.xx(1)-0.5.*lx; y0=popCont.yy(1)-0.5.*ly;
            
            indx=ceil((xxJF-x0)./lx); indy=ceil((yyJF-y0)./ly);
            outScore=ones(size(indx));
            for ii=1:length(xxJF)
                outScore(ii)=popCont.popContour(indy(ii),indx(ii));
            end

            %% plot 
            
            figure('position',[1432 421 1000 750],'Color','w')
            
            %% underlying hist 

            hh=scatterhist(xx,yy,'Group',jfscore,...
                'location','northeast','direction','out','legend','off','color',cols,...
                'markersize',1);
            
            %0.6950    0.6960
            set(hh(1),'position',[0.1000    0.1000   0.7 0.7 ],'fontsize',14);
            hh(1).YLabel.String='';
            hh(1).XLabel.String='';
            %'XTick',[],'YTick',[],'YAxisLocation','left','XAxisLocation','bottom');
            set(hh(2),'position',[0.80 0.1 0.119806451612903 0.75]);
            set(hh(3),'position',[0.1 0.80 0.75 0.132374233128834]);
            hh(2).Children(1).LineWidth=1.5; hh(2).Children(2).LineWidth=1.5;
            hh(3).Children(1).LineWidth=1.5; hh(3).Children(2).LineWidth=1.5;
            xlim(xl);ylim(yl);
            
            %% map 
            
            ax1=axes;
            axPos=get(hh(1),'position');
            set(ax1,'position',axPos);
            
            imagesc(xl,yl,squeeze(bird(:,:,1)))
            set(gca,'ydir','normal','fontsize',14)
            
            colormap(cmap)
            
            %% JF contour 
            
            hold on
            [~,h(1)]=contour(popCont.xx,popCont.yy,popCont.popContour,'ShowText','off','LineColor',[0 0 1],...
                'LineColor',cols(2,:),'linewidth',2,...
                'LevelList',[98 75:-25:5],'Fill','off','linestyle','-',...
                'DisplayName','Jellyfish');
            
            
            plot(xxJF(outScore>98),yyJF(outScore>98),'^',...
                'color',cols(2,:),'markersize',6.5,...
                'linewidth',1.5,'MarkerFaceColor',otherCol(2,:),...
                'Displayname','JF beyond 98 percentile');
            grid
            %legend(h)
            
            xfac=0.83; yfac=0.08;
            text(xfac.*diff(xl)+xl(1),yfac.*diff(yl)+yl(1),sims{k},...%'Edgecolor','k','backgroundcolor',[1,0.97,0.97],...
                'Interpreter','latex','fontsize',17,'fontweight','bold','color','k')
            
            xlim(xl);ylim(yl);

            
            xlabelmine(xlab{i},16);
            ylabelmine(ylab{j},16);


            linkaxes([hh(1),ax1]) 
            
            %% print figure 
            fname=sprintf('jfProps_%s_%s_%s',xfields{i},yfields{j},sims{k});
            if  printFlag; printout_fig(gcf,fname,'nopdf','v'); end
    
            
        end
        close all
    end
    
    
end