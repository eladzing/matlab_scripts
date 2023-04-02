%% generate 5 score ranges

ranges=0:0.1:1.01;
ranges(end)=ranges(end)+0.1;
global DEFAULT_PRINTOUT_DIR
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/mosaics'];

nn=10;

%% get 20 objects of a given score range
for i=1:length(ranges)-1
    
    inds=find(objectTable.scoreWeighted>=ranges(i) & ...
        objectTable.scoreWeighted<ranges(i+1)) ;
    
    % select nn objects.
    inds=shuffleArray(inds);
    inds=inds(1:nn);
    imageStruct(i).tags=objectTable.tag(inds);
    imageStruct(i).score=objectTable.scoreWeighted(inds);
end

%% parse tags and read in images.

fnameBase='%s_subhalo%s_snap%s_map_conf0_contours_conf1.png';
baseImagepath='\\wsl.localhost\Ubuntu-20.04\home\kipod\sshfsMounts\vera\IllustrisTNG\simulationOutput\';
for k=1:length(imageStruct)
    
    for i=1:length(imageStruct(k).tags)
        
        if imageStruct(k).tags(i).contains("TNG50")
            simTag='L35n2160TNG';
        elseif imageStruct(k).tags(i).contains("TNG100")
            simTag='L75n1820TNG';
        else
            error('what?')
        end
        
        imagePath=[baseImagepath simTag '\postprocessing\Zooniverse_CosmologicalJellyfish\images_all\'];
        
        snap=imageStruct(k).tags(i).extractBetween('snp','subid');
        
        snapPath=join([imagePath simTag '_snap' snap '\'],'');
        
        subid=imageStruct(k).tags(i).extractBetween('subid','typ');
        fname=sprintf(fnameBase,simTag,subid,snap);
        fullname=join([snapPath fname],'');
        imageStruct(k).image{i}=imread(fullname);
    end
end
%% plot 5 X 4 image page

rangesTag=({'0--0.2','0.2--0.4','0.4--0.6','0.6--0.8','0.8--1.0'});

for k=1:2:length(imageStruct)
    
    hf=myFigure('pos',[ 1425 70  1000 1290]);
    
    tt=tiledlayout(5,4);
    
    for i=1:length(imageStruct(k).image)
        
        nexttile
        imagesc(imageStruct(k).image{i});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
        scoreTag=sprintf('score=%3.2f',imageStruct(k).score(i));
        text(100,200,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
        
    end
    
    for i=1:length(imageStruct(k+1).image)
        
        nexttile
        imagesc(imageStruct(k+1).image{i});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
        scoreTag=sprintf('score=%3.2f',imageStruct(k+1).score(i));
        text(100,200,scoreTag,'fontsize',16,'color','w','Interpreter','latex');
        
    end
    
    
    
    set(tt,'Padding','none','TileSpacing','none')
    tt.Title.String=['score range ' rangesTag{ceil(k/2)}];
    tt.Title.Interpreter='latex';
    tt.Title.FontSize=20;
    
    
    outname=sprintf('cjf_mosaic_appendix_scoreBin_%i',ceil(k/2));
    if printFlag; printout_fig(gcf,outname,'nopdf','v','printoutdir',outdir); end
    
end


