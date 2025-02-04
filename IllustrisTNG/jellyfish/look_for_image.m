
global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '\cjf_image_mosiac_2.mat'])
global DEFAULT_PRINTOUT_DIR
outdir=[DEFAULT_PRINTOUT_DIR '/jellyfish/paper'];

%% set indices
inds=[1 7;1 5;...
    3 5; 4 9;...
    5 8; 7 2;...
    9 7; 7 1;...
    3 8; 5 10];


%% plot 5 X 4 image page
ranges=fliplr({'0--0.2','0.2--0.4','0.4--0.6','0.6--0.8','0.8--1.0'});

%hf=myFigure('pos',[ 1425 70  1050 1250]);

%tt=tiledlayout(5,4);
%cnt=1;



%cnt2=1;

for i=1:length(imageStruct)
    
    for j=1:length(imageStruct(i).score)
        %ii=inds(i,j);
        %nexttile
        figure(gcf)
        imagesc(imageStruct(i).image{j});
        axis square
        set(gca,'Ytick',[],'xtick',[],'box','off');
        
        scoreTag=sprintf('score=%3.2f',imageStruct(i).score(j));
        text(100,150,scoreTag,'fontsize',40,'color','w','Interpreter','latex');
        %if mod(cnt,4)==1
        %    ylabelmine(ranges{cnt2},20);
        %        cnt2=cnt2+1;
        %end
    %cnt=cnt+1;
    fprintf([ num2str(i) '/' num2str(j)  ' ']);
    %pause
    end
    
end


% set(tt,'Padding','none','TileSpacing','none')
% tt.YLabel.String='Score Range';
% tt.YLabel.Interpreter='latex';
% tt.YLabel.FontSize=18;
%% print
outname=sprintf('cjf_scoreMosaic');
if printFlag; printout_fig(gcf,outname,'pdf','v','printoutdir',outdir); end


