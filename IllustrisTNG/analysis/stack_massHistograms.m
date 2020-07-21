
%% load stuff
global DEFAULT_MATFILE_DIR
global simDisplayName
global DRACOFLAG
snap=99;

if readFlag
    fname=sprintf('gasProperties_massHistograms_snp%s_%s.mat',num2str(snap),simDisplayName);
    load([DEFAULT_MATFILE_DIR '/' fname],'massHist')
    
    if DRACOFLAG
        fofs=illustris.groupcat.loadHalos(bp,snap);
        subs=illustris.groupcat.loadSubhalos(bp,snap);
    else
        loadFofSubTNG100
    end
end

%% arrange mask
if maskFlag
    subsInfo=illustris.infrastructure.build_sub_fof_connection(subs,fofs);
    spbMask = mk_splashback_mask('time',5,'both',0.1);
    massThresh=10^9;
    galMask=illustris.utils.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'snap',snap);

    totMask= subsInfo.isCentral  & ~spbMask;
    totMask=totMask(massHist.galMask);
    
    %massRange=9:0.5:12.5;
    massRange=9:12;
    
    galMass=massHist.galMass(massHist.galMask);
    galMass=galMass(totMask);
    mInd=discretize(log10(galMass),massRange);
    
    indxs=find(totMask);
    
    ssfr=illustris.utils.calc_ssfr(subs);
    ssfr=ssfr(massHist.galMask);
    qMask0=ssfr<1e-11;
    
    ssfr=ssfr(totMask);
    qMask=ssfr<1e-11;
end
%% set framework for global histograms


if buildFlag
    nb=30;
    
    hType={'dens','ent','temp','tcTff','tc'};
    
    
    fld={'CGM','Out','Sub'};%,'Gal'
    
    fieldRange.tc=linspace(-0.5,1.5,nb+1);
    fieldRange.tcTff=linspace(-2,2,nb+1);
    fieldRange.temp=linspace(4,8,nb+1);
    fieldRange.ent=linspace(-3,3,nb+1);
    fieldRange.dens=linspace(-5,0,nb+1);
    
    for i=1:length(hType)
        histos.Gal.(hType{i})=zeros(length(massRange)-1,nb+2);
        histos.Gal.([hType{i} 'Q1'])=zeros(length(massRange)-1,nb+2);
        histos.Gal.([hType{i} 'SF1'])=zeros(length(massRange)-1,nb+2);
        histos.Gal.([hType{i} 'Q2'])=zeros(length(massRange)-1,nb+2);
        histos.Gal.([hType{i} 'SF2'])=zeros(length(massRange)-1,nb+2);
    end
    
    
    histos.CGM=histos.Gal;
    histos.Out=histos.Gal;
    histos.Sub=histos.Gal;
    
    %% run over mass ranges
    nSampleBase=100;
    
    for im=1:length(massRange)-1
        
        inds=indxs(mInd==im);
        qm=qMask(mInd==im);
        
        indsQ=inds(qm);
        indsS=inds(~qm);
        
        %scramble list
        inds=shuffleArray(inds);
        indsQ=shuffleArray(indsQ);
        indsS=shuffleArray(indsS);
        
        nSample(im)=min(nSampleBase,length(inds));
        nSampleQ(im)=min(nSampleBase,length(indsQ));
        nSampleS(im)=min(nSampleBase,length(indsS));
        
        
        % run over different regions
        for k=1:length(fld)
            
            % run over different parameters
            for j=2 %  1:length(hType)
                
                % run over sample (total and sub-divide to Q and SF)
                cntQ=0;
                cntS=0;
                for i=1:nSample(im)
                    
                    ind=inds(i);
                    
                    % get mass hisogram
                    xx=massHist.(['in' fld{k}]).(hType{j})(:,ind);
                    yy=massHist.(['in' fld{k}]).([hType{j} 'Mass'])(:,ind);
                    if sum(yy)>0
                        yy=yy./sum(yy); % fractional mass - good for stacking
                    end
                    
                    % distribute
                    hInd=discretize(xx,fieldRange.(hType{j}))+1;
                    hInd(xx<fieldRange.(hType{j})(1))=1;
                    hInd(xx>fieldRange.(hType{j})(end))=nb+2;
                    
                    for ii=1:nb+2
                        histos.(fld{k}).(hType{j})(im,ii)=histos.(fld{k}).(hType{j})(im,ii)+...
                            sum(yy(hInd==ii));
                    end
                    
                    if qMask0(ind)
                        cntQ=cntQ+1;
                        
                        for ii=1:nb+2
                            histos.(fld{k}).([hType{j} 'Q1'])(im,ii)=histos.(fld{k}).([hType{j} 'Q1'])(im,ii)+...
                                sum(yy(hInd==ii));
                        end
                    else
                        cntS=cntS+1;
                        for ii=1:nb+2
                            histos.(fld{k}).([hType{j} 'SF1'])(im,ii)=histos.(fld{k}).([hType{j} 'SF1'])(im,ii)+...
                                sum(yy(hInd==ii));
                        end
                    end
                    
                end
                %bar(1:nb+2,histos.(fld{k}).(hType{j})(im,:))
                
                
                
                histos.(fld{k}).(hType{j})(im,:)=histos.(fld{k}).(hType{j})(im,:)./nSample(im);
                if cntS>0
                    histos.(fld{k}).([hType{j} 'SF1'])(im,:)= histos.(fld{k}).([hType{j} 'SF1'])(im,:)./cntS;
                end
                if cntQ>0
                    histos.(fld{k}).([hType{j} 'Q1'])(im,:)= histos.(fld{k}).([hType{j} 'Q1'])(im,:)./cntQ;
                end
                
                nSampleQ1(im)=cntQ;
                nSampleS1(im)=cntS;
                
                % run over Quenched
                for i=1:nSampleQ(im)
                    
                    ind=indsQ(i);
                    
                    % get mass hisogram
                    xx=massHist.(['in' fld{k}]).(hType{j})(:,ind);
                    yy=massHist.(['in' fld{k}]).([hType{j} 'Mass'])(:,ind);
                    if sum(yy)>0
                        yy=yy./sum(yy); % fractional mass - good for stacking
                    end
                    
                    if strcmp(hType{j},'ent')
                        plot(xx,yy)
                        titlemine(sprintf('m%i %s Q, %i out of %i',...
                            massRange(im),fld{k},i,nSampleQ(im)));
                        pause
                    end
                    
                    
                    
                    % distribute
                    hInd=discretize(xx,fieldRange.(hType{j}))+1;
                    hInd(xx<fieldRange.(hType{j})(1))=1;
                    hInd(xx>fieldRange.(hType{j})(end))=nb+2;
                    
                    for ii=1:nb+2
                        histos.(fld{k}).([hType{j} 'Q2'])(im,ii)=histos.(fld{k}).([hType{j} 'Q2'])(im,ii)+...
                            sum(yy(hInd==ii));
                    end
                    
                end
                
                % run over Star-Forming
                for i=1:nSampleS(im)
                    
                    ind=indsS(i);
                    
                    % get mass hisogram
                    xx=massHist.(['in' fld{k}]).(hType{j})(:,ind);
                    yy=massHist.(['in' fld{k}]).([hType{j} 'Mass'])(:,ind);
                    if sum(yy)>0
                        yy=yy./sum(yy); % fractional mass - good for stacking
                    end
                    
                      if strcmp(hType{j},'ent')
                        plot(xx,yy)
                        titlemine(sprintf('m%i %s SF, %i out of %i',...
                            massRange(im),fld{k},i,nSampleQ(im)));
                        pause
                    end
                    
                    % distribute
                    hInd=discretize(xx,fieldRange.(hType{j}))+1;
                    hInd(xx<fieldRange.(hType{j})(1))=1;
                    hInd(xx>fieldRange.(hType{j})(end))=nb+2;
                    
                    for ii=1:nb+2
                        histos.(fld{k}).([hType{j} 'SF2'])(im,ii)=histos.(fld{k}).([hType{j} 'SF2'])(im,ii)+...
                            sum(yy(hInd==ii));
                    end
                    
                end
                
                
                
                
                
                histos.(fld{k}).([hType{j} 'SF2'])(im,:)=...
                    histos.(fld{k}).([hType{j} 'SF2'])(im,:)./nSampleS(im);
                histos.(fld{k}).([hType{j} 'Q2'])(im,:)=...
                    histos.(fld{k}).([hType{j} 'Q2'])(im,:)./nSampleQ(im);
                
            end
        end
    end
end



%% plot
if plotFlag
    
    col1=brewermap(8,'Set1');
    col2=brewermap(8,'Set2');
    col3=brewermap(8,'Paired');
    %nn=1:nb+2;
    hLab={'log Density $[\mathrm{cm^{-2}}]$',...
        'Entropy $[\mathrm{Kev\, cm^{2}}]$',...
        'log Temperature $[\mathrm{K}]$',...
        'log $t_\mathrm{cool}/t_\mathrm{ff}$',...
        'log Cooling Time $[\mathrm{Gyr}]$'};
    
    fldTag={'Inner CGM' 'Outer CGM' 'Total CGM' 'Galactic Gas'};
    %msLab={'9_95' '95_10' '10_105' '105_11' '11_115' '115_12' '12'};
    msLab={'9_10' '10_11' '11_12'};
    
    for im=1:length(massRange)-1
        
        %         mStr=sprintf('$\\mathrm{M_{star}}=[%s,%s] \\,\\mathrm{M_\\odot}$',...
        %             num2str(massRange(im)),num2str(massRange(im+1)));
        mStr=sprintf('$[%s,%s] \\,\\mathrm{M_\\odot}$',...
            num2str(massRange(im)),num2str(massRange(im+1)));
        % run over different regions
        for k=1:length(fld)
            
            % run over different parameters
            for j=1:length(hType)
                clear xxx
                
                xxx(1:nb+1)=fieldRange.(hType{j});
                xxx(nb+2)=xxx(nb+1)+diff(xxx(nb:nb+1));
                lend=1:nb+2;
                
                %                 if j~=4
                %                     lend=1:nb+2;
                %                 else
                %                     %    xxx(1:nb+1)=fieldRange.(hType{j});
                %                     lend=2:nb+1;
                %                 end
                dx=0.25.*max(diff(xxx));
                
                bb=100.*histos.(fld{k}).(hType{j})(im,lend);
                
                %             bbQ=100.*histos.(fld{k}).([hType{j} 'Q1'])(im,1:lend);
                %             bbS=100.*histos.(fld{k}).([hType{j} 'SF1'])(im,1:lend);
                %
                %             bbb=cat(2,bb',bbQ',bbS');
                %             bbb(isnan(bbb))=0;
                %             %% Q/SF subdivision
                %             figure('color','w')
                %
                %             b=bar(xxx,bbb);
                %
                %             b(1).FaceColor='k';
                %             b(2).FaceColor=cols(1,:);
                %             b(3).FaceColor=cols(2,:);
                %
                %             b(1).DisplayName=['All (' num2str(nSample(im)) ')'];
                %             b(2).DisplayName=['Q (' num2str(nSampleQ1(im)) ')'];
                %             b(3).DisplayName=['SF (' num2str(nSampleS1(im)) ')'];
                %
                %             hl=legend;
                %             set(hl,'Interpreter','latex','fontsize',14);
                %
                %             xlabelmine(hLab{j},16);
                %             ylabelmine('percent',16);
                %
                %
                %             titlemine(sprintf('%s, %s',fld{k},mStr));
                %
                %             set(gca,'fontsize',14)
                %
                %             fname=sprintf('gasHistogramSub_%s_%s_mass%s_nSamp%s_%s',...
                %                 fld{k},hType{j},msLab{im},num2str(nSample(im)),simDisplayName);
                %           %  printout_fig(gcf,fname,'subdir','gasHistograms','v')
                %
                %             close(gcf)
                
                %% Q/SF sample
                
                bbQ=100.*histos.(fld{k}).([hType{j} 'Q2'])(im,lend);
                bbS=100.*histos.(fld{k}).([hType{j} 'SF2'])(im,lend);
                
                bbb=cat(2,bb',bbQ',bbS');
                bbb(isnan(bbb))=0;
                
                figure('color','w')
                
                
                %             h(1)=stairs(xxx,bbb(:,1),'linewidth',1.5,'color','k',...
                %                'DisplayName',['All (' num2str(nSample(im)) ')']);
                %
                
                h(1)=stairs(xxx(lend),bbb(:,2),'linewidth',1.8,'color',col1(1,:),...
                    'DisplayName',['Q (' num2str(nSampleQ(im)) ')']);
                hold on
                h(2)=stairs(xxx(lend)-dx,bbb(:,3),'linewidth',1.8,'color',col3(2,:),...
                    'DisplayName',['SF (' num2str(nSampleS(im)) ')']);
                
                xlim([xxx(1) xxx(end)])
                
                hl=legend(h);
                set(hl,'Interpreter','latex','fontsize',14);
                
                xlabelmine(hLab{j});
                ylabelmine('percent');
                
                set(gca,'fontsize',14)
                %titlemine(sprintf('%s, %s',fld{k},mStr));
                xl=xlim;
                yl=ylim;
                
                text(xl(1)+0.05*diff(xl),yl(2)-0.05*diff(yl),mStr,...
                    'Interpreter','latex','Fontsize',14)
                text(xl(1)+0.3*diff(xl),yl(2)-0.05*diff(yl),fldTag{k},...
                    'Interpreter','latex','Fontsize',14)
                
                
                
                
                
                fname=sprintf('gasHistogramQSF_%s_%s_mass%s_nSamp%s_%s',...
                    fld{k},hType{j},msLab{im},num2str(nSample(im)),simDisplayName);
                printout_fig(gcf,fname,'subdir','gasHistograms','v')
                
                close(gcf)
                
            end
        end
    end
    
end



%             b=bar(xxx,bbb);
%
%             b(1).FaceColor='k';
%             b(2).FaceColor=cols(1,:);
%             b(3).FaceColor=cols(2,:);
%
%             b(1).DisplayName=['All (' num2str(nSample(im)) ')'];
%             b(2).DisplayName=['Q (' num2str(nSampleQ(im)) ')'];
%             b(3).DisplayName=['SF (' num2str(nSampleS(im)) ')'];
%             b(1).BarWidth=1;
%             b(2).BarWidth=1;
%             b(3).BarWidth=1;
%             hl=legend;









