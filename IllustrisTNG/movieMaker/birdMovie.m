function makeBirdMovie(galhist)


hf=figure;

set(hf,'position',[220 605 1642 718]);




cnt=1;

for i=75:-1:1
    


bird1=squeeze(galhist.inGal.bird(:,:,i));
bird2=squeeze(galhist.inSub.bird(:,:,i));

subplot(1,2,1)
illustris.plots.plot_phaseDiagram(galhist.birdXlim,galhist.birdYlim,...
    bird1,'fig',hf,'caxis',[7.5 10]);
titlemine(sprintf('Gal $z=%4.2f$',galhist.zred(i)));

subplot(1,2,2)
illustris.plots.plot_phaseDiagram(galhist.birdXlim,galhist.birdYlim,...
    bird2,'fig',hf,'caxis',[7.5 10]);
titlemine(sprintf('Sub $z=%4.2f$',galhist.zred(i)));
%fprintf('i= %i \n',i)


F(cnt)=getframe(gcf);



cnt=cnt+1;


end

global DEFAULT_MOVIE_DIR

Name=

vid = VideoWriter('bird3');
set(vid,'FrameRate',5);
open(vid)

for i=1:length(F)
    
  writeVideo(vid,F(i));
end


 close(vid)
 
 
 
 % set properties 
 
 
 
 

 
 
 
 
% 
% writeVideo(vid,F)