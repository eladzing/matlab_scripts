%% look for objects which afe JF in optimized and non-Jf in random




global DEFAULT_PRINTOUT_DIR
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/mosaics'];

fnameBaseR='%s_subhalo%s_snap%s_map_conf0_contours_conf1.png';
fnameBaseP='%s_subhalo%s_snap%s_map_conf0_contours_conf1_SubhaloVel_in_image_plane.png';
baseImagepath='\\wsl.localhost\Ubuntu-20.04\home\kipod\sshfsMounts\vera\IllustrisTNG\simulationOutput\';

%% get  objects of a given score range
inds=find(objectTableComp.scoreWeightedPref>=0.8 & ...
    objectTableComp.scoreWeightedRand<0.7);

for i=1:length(inds)
    ii=inds(i);
    tag=objectTableComp.tag(ii);
    imageStruct.tags(i)=tag;
    imageStruct.scoreRand(i)=objectTableComp.scoreWeightedRand(ii);
    imageStruct.scorePref(i)=objectTableComp.scoreWeightedPref(ii);
    
    if tag.contains("TNG50")
        simTag='L35n2160TNG';
    elseif tag.contains("TNG100")
        simTag='L75n1820TNG';
    else
        error('what?')
    end
    
    snap=tag.extractBetween('snp','subid');
    subid=tag.extractBetween('subid','typ');
    
    imagePath=[baseImagepath simTag '\postprocessing\Zooniverse_CosmologicalJellyfish\images_all\'];
    
    
    snapPathRand=join([imagePath simTag '_snap' snap '\'],'');
    snapPathPref=join([imagePath simTag '_snap' snap '_SubhaloVel_in_image_plane\'],'');
    
    fnameRand=sprintf(fnameBaseR,simTag,subid,snap);
    fnamePref=sprintf(fnameBaseP,simTag,subid,snap);
    fullnameRand=join([snapPathRand fnameRand],'');
    fullnamePref=join([snapPathPref fnamePref],'');
    imageStruct.imageRand{i}=imread(fullnameRand);
    imageStruct.imagePref{i}=imread(fullnamePref);
    
    
    
end

%% plot 5 X 4 image page
cnt=1;
figCnt=1;
while cnt<=length(inds)
    
    hf=myFigure('pos',[ 1425 70  1000 1250]);
    
    tt=tiledlayout(5,4);
    
    for i=1:10
        
        if cnt>length(inds)
            break
        end
        
        nexttile
        imagesc(imageStruct.imageRand{cnt});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
        scoreTag=sprintf('%3.2f R',imageStruct.scoreRand(cnt));
        text(100,200,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
        text(100,350,imageStruct.tags(cnt).extractBefore('typ'),'fontsize',11,'color','w');

        nexttile
        imagesc(imageStruct.imagePref{cnt});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
        
        %idTag=sprintf('%s ',imageStruct.tags(cnt));
        scoreTag=sprintf('%3.2f P',imageStruct.scorePref(cnt));
        text(100,200,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
        text(100,350,imageStruct.tags(cnt).extractBefore('typ'),'fontsize',11,'color','w');
        cnt=cnt+1;
    end
    
    
    set(tt,'Padding','none','TileSpacing','none')
    
    outname=sprintf('cjf_mosaic_randNon_optJF_%i',figCnt);
    if printFlag; printout_fig(gcf,outname,'nopdf','v','printoutdir',outdir); end
    figCnt=figCnt+1;
end


%% Repeat for objects which are JF in random and non-jf in optimized
inds=find(objectTableComp.scoreWeightedPref<0.8 & ...
    objectTableComp.scoreWeightedRand>=0.8);

for i=1:length(inds)
    ii=inds(i);
    tag=objectTableComp.tag(ii);
    imageStruct2.tags(i)=tag;
    imageStruct2.scoreRand(i)=objectTableComp.scoreWeightedRand(ii);
    imageStruct2.scorePref(i)=objectTableComp.scoreWeightedPref(ii);
    
    if tag.contains("TNG50")
        simTag='L35n2160TNG';
    elseif tag.contains("TNG100")
        simTag='L75n1820TNG';
    else
        error('what?')
    end
    
    snap=tag.extractBetween('snp','subid');
    subid=tag.extractBetween('subid','typ');
    
    imagePath=[baseImagepath simTag '\postprocessing\Zooniverse_CosmologicalJellyfish\images_all\'];
    
    
    snapPathRand=join([imagePath simTag '_snap' snap '\'],'');
    snapPathPref=join([imagePath simTag '_snap' snap '_SubhaloVel_in_image_plane\'],'');
    
    fnameRand=sprintf(fnameBaseR,simTag,subid,snap);
    fnamePref=sprintf(fnameBaseP,simTag,subid,snap);
    fullnameRand=join([snapPathRand fnameRand],'');
    fullnamePref=join([snapPathPref fnamePref],'');
    imageStruct2.imageRand{i}=imread(fullnameRand);
    imageStruct2.imagePref{i}=imread(fullnamePref);
end

%% plot 5 X 4 image page
cnt=1;
figCnt=1;
while cnt<=length(inds)
    
    hf=myFigure('pos',[ 1425 70  1000 1250]);
    
    tt=tiledlayout(5,4);
    
    for i=1:10
        
        if cnt>length(inds)
            break
        end
        
        nexttile
        imagesc(imageStruct2.imageRand{cnt});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
         scoreTag=sprintf('%3.2f R',imageStruct2.scoreRand(cnt));
        text(100,200,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
        text(100,350,imageStruct2.tags(cnt).extractBefore('typ'),'fontsize',11,'color','w');
        
        nexttile
        imagesc(imageStruct2.imagePref{cnt});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
        scoreTag=sprintf('%3.2f P',imageStruct2.scorePref(cnt));
        text(100,200,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
        text(100,350,imageStruct2.tags(cnt).extractBefore('typ'),'fontsize',11,'color','w');
        
        cnt=cnt+1;
    end
    
    
    set(tt,'Padding','none','TileSpacing','none')
    
    outname=sprintf('cjf_mosaic_randJF_optNon_%i',figCnt);
    if printFlag; printout_fig(gcf,outname,'nopdf','v','printoutdir',outdir); end
    figCnt=figCnt+1;
end