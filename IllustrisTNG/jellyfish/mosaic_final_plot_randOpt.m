
global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '\cjf_image_mosiac_randOpt.mat'])
global DEFAULT_PRINTOUT_DIR
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/mosaics'];




%%  opt JF / rand non : set indices
list1=["TNG50snp067subid00000098" ...
    "TNG50snp067subid00115166" ...
    "TNG50snp067subid00139885" ...
    "TNG50snp067subid00202863" ...
    "TNG50snp099subid00000037" ...
    "TNG50snp099subid00220620" ...
    "TNG50snp099subid00429476" ...
    "TNG100snp067subid00020427" ...
    "TNG100snp067subid00007390" ...
    "TNG100snp067subid00166403"];

list2=["TNG50snp099subid00229944" ...
    "TNG50snp099subid00313699" ...
    "TNG100snp067subid00152090" ...
    "TNG100snp067subid00020427"];
list=list2;

list3=["TNG50snp067subid00115161" ...
    "TNG100snp099subid00306062" ...
    "TNG100snp099subid00242223" ...
    "TNG50snp099subid00096797"...
    "TNG100snp067subid00069163"...
    "TNG100snp067subid00078691"...
    "TNG100snp099subid00206230"];
list=list3;

list4=["TNG100snp099subid00381151" ...
    "TNG100snp099subid00083332" ...
    "TNG50snp067subid00077304" ...
    "TNG50snp067subid00115165" ...
    "TNG50snp067subid00139881" ...
    "TNG50snp067subid00258172"];


%list=list;



%% plot 2 X 4 image page
%ranges=fliplr({'0--0.2','0.2--0.4','0.4--0.6','0.6--0.8','0.8--1.0'});

%%
list=list1;
sameTag='_same';

imsgtr=imageStruct;

ii=find(imsgtr.tags.contains(list));
tag=['randN_OptJF' sameTag];
if isempty(ii)
    imsgtr=imageStruct2;
    tag=['randJF_OptN' sameTag];
    ii=find(imsgtr.tags.contains(list));
    if isempty(ii)
        error('%s - wrong list or structure',current_function().upper);
    end
end



for i=1:length(list)
    
    hf=myFigure('pos',[ 666  325 1800 900]);
    
    tt=tiledlayout(1,2);
    
    ii=find(imsgtr.tags.contains(list(i)));
    nexttile
    imagesc(imsgtr.imageRand{ii});
    axis square
    set(gca,'Ytick',[],'xtick',[],'box','off');
    
    scoreTag=sprintf('score=%3.2f',imsgtr.scoreRand(ii));
    text(100,200,'Random','fontsize',32,'color','w','Interpreter','latex');
    text(100,350,scoreTag,'fontsize',32,'color','w','Interpreter','latex');
    
    nexttile
    imagesc(imsgtr.imagePref{ii});
    axis square
    set(gca,'Ytick',[],'xtick',[],'box','off');
    
    scoreTag=sprintf('score=%3.2f',imsgtr.scorePref(ii));
    text(100,200,'Optimized','fontsize',32,'color','w','Interpreter','latex');
    text(100,350,scoreTag,'fontsize',32,'color','w','Interpreter','latex');
    
    
    
    
    set(tt,'Padding','none','TileSpacing','none')
    % tt.YLabel.String='Score Range';
    % tt.YLabel.Interpreter='latex';
    % tt.YLabel.FontSize=18;
    %% print
    outname=sprintf('cjf_%s_mosaic_%i',tag,i);
    %outname=sprintf('cjf_randOpt_same2_mosaic_%i',i);
    if printFlag; printout_fig(gcf,outname,'pdf','v','printoutdir',outdir); end
    pause
end
%%  opt non / rand JFn : set indices

%
%
%
% %% plot 2 X 4 image page
% %ranges=fliplr({'0--0.2','0.2--0.4','0.4--0.6','0.6--0.8','0.8--1.0'});
%
% hf=myFigure('pos',[ 666  325 1200 600]);
%
% tt=tiledlayout(2,4);
%
% for i=1:length(inds)
%
%
%     ii=inds(i);
%     nexttile
%     imagesc(imageStruct2.imageRand{ii});
%     axis square
%     set(gca,'Ytick',[],'xtick',[],'box','off');
%
%     scoreTag=sprintf('score=%3.2f',imageStruct2.scoreRand(ii));
%     text(100,200,'Random','fontsize',16,'color','w','Interpreter','latex');
%     text(100,350,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
%
%     nexttile
%     imagesc(imageStruct2.imagePref{ii});
%     axis square
%     set(gca,'Ytick',[],'xtick',[],'box','off');
%
%     scoreTag=sprintf('score=%3.2f',imageStruct2.scorePref(ii));
%     text(100,200,'Optimized','fontsize',16,'color','w','Interpreter','latex');
%     text(100,350,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
%
% end
%
%
% set(tt,'Padding','none','TileSpacing','none')
% % tt.YLabel.String='Score Range';
% % tt.YLabel.Interpreter='latex';
% % tt.YLabel.FontSize=18;
% %% print
% outname=sprintf('cjf_randOpt_Mosaic2');
% if printFlag; printout_fig(gcf,outname,'nopdf','v','printoutdir',outdir); end
%
%
%
