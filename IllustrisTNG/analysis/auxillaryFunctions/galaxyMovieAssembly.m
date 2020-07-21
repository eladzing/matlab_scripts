function galaxyMovieAssembly(F,name,varargin)
%GALAXYMOVIEASSEMBLY Auxilary function for converting frames into movie

frameRate=2;
mpegFlag=true;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'frame','frameRate','rate','fr'}
            i=i+1;
            frameRate=varargin{i};
        case {'avi','nompfeg'}
            mpegFlag=false;
        otherwise
            str=upper(galaxyMovieAssembly);
            error('%s - Illegal argument: %s',str,varargin{i});
    end
    i=i+1;
end



if ~exist('frameRate','var')
    frameRate=2;
end

global DEFAULT_MOVIE_DIR

movieName=[DEFAULT_MOVIE_DIR '/' name];

vid = VideoWriter(movieName);

set(vid,'FrameRate',frameRate);
open(vid)

for i=1:length(F)
    
    writeVideo(vid,F(i));
end

close(vid)

%% convert to mp4
if mpegFlag
    %comm=sprintf('ffmpeg -i %s.avi -acodec libfaac -b:a 128k -vcodec mpeg4 -b:v 1200k -flags +aic+mv4 %s.mp4',movieName,movieName);
        comm=sprintf('ffmpeg -i %s.avi %s.mp4',movieName,movieName);
        
    system(comm);
    system(sprintf('rm -f %s.avi',movieName));
    
end

