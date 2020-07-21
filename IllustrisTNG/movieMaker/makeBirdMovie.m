function makeBirdMovie(galhist,name)


%% prepare movie frames 

hf=figure;

set(hf,'position',[220 605 1642 718]);


cnt=1;

for i=length(galhist.zred):-1:1
    


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

%% write movie to file 

global DEFAULT_MOVIE_DIR

movieName=[DEFAULT_MOVIE_DIR '/' name];

vid = VideoWriter(movieName);

set(vid,'FrameRate',5);
open(vid)

for i=1:length(F)
    
  writeVideo(vid,F(i));
end


 close(vid)
 
 %% convert to mp4 
 
 comm=sprintf('ffmpeg -i %s.avi -acodec libfaac -b:a 128k -vcodec mpeg4 -b:v 1200k -flags +aic+mv4 %s.mp4',movieName,movieName);
 
 system(comm);
 system(sprintf('rm -f %s.avi',movieName));
 
 
 
 
 
 % set properties 
 
 
 
 

 
 
 
 
% 
% writeVideo(vid,F)