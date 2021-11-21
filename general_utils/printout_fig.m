function printout_fig(figHandle,name,varargin)

%printoutDirBase='C:\Users\eladzing\Documents\cluster\printout';

global DEFAULT_PRINTOUT_DIR

printoutDir=DEFAULT_PRINTOUT_DIR;
%format='both';
%format='pdf';

verboseflag=false ;
figFlag=false;
epsFlag=false;
pdfFlag=false;
pngFlag=true;

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'verbose','v'}
            verboseflag=true;
        case {'dir','printout','printoutdir'}
            i=i+1;
            printoutDir=varargin{i};
        case {'subdir'}
            i=i+1;
            dir=varargin{i};
            printoutDir=sprintf('%s/%s',printoutDir,dir);
        case {'format','type'}
            
            error('printout_fig: format is obsolete. please fix script accordingly');
            %format=varargin{i};  OBSOLETE
        case {'png'}
            pngFlag=true;
        case {'eps'}
            epsFlag=true;
        case {'pdf'}
            pdfFlag=true;
        case {'nopdf'}
            pdfFlag=false;
        case {'nopng'}
            pngFlag=false;
        case {'nofig'}
            figFlag=false;
        case {'fig','mfig'}
            figFlag=true;
        otherwise
            error('printout_fig:illegal option:%s',varargin{i} );
    end
    i=i+1;
end

%  switch(format)
%      case('png')
%          printformat=[true false];
%      case('eps')
%          printformat=[false true];
%      case('all')
%          printformat=[true true];
%      otherwise
% %         error('printout_fig:illegal format');
% % end

%test
% set(figHandle,'PaperUnits', 'inches');
%
% width = 5;     % Width in inches
% height = 5;    % Height in inches
% alw = 0.75;    % AxesLineWidth
% fsz = 12;      % Fontsize
% lw = 2;      % LineWidth
% msz = 8;       % MarkerSize
%
% pos = get(figHandle, 'Position');
% set(figHandle, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

% set(figHandle,'InvertHardcopy','on');
%
% papersize = get(figHandle, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(figHandle,'PaperPosition', myfiguresize);

%% check to see if printout directory exists exists
dirStat=exist(printoutDir,'dir');
if dirStat==0
    mkdir(printoutDir)
    dirStat=exist(printoutDir,'dir');
    
    if dirStat~=7
        error('PRINTOUT_FIG - something wrong with directory: %s',printoutDir);
    end
end
   
%% 

figure(figHandle)
fulln=sprintf('%s/%s',printoutDir,name);

% PNG version
if pngFlag
    fullname=sprintf('%s/%s.%s',printoutDir,name,'png');
    if verboseflag
        fprintf('printing to %s \n',fullname);
    end
    %exportfig(figHandle,fullname,'format','png','color','rgb','width',20);%,'Resolution',300,'Fontmode','scaled');
    export_fig(fullname,'-png')
end

%figure(figHandle);
%tightfig(figHandle);
% EPS version
if epsFlag
    fullname=sprintf('%s/%s.%s',printoutDir,name,'eps');
    %fulln=sprintf('%s/%s',printoutDir,name);
    if verboseflag
        fprintf('printing to %s \n',fullname);
    end
    
    export_fig(fullname,'-eps','-transparent')
end

% PDF version
if pdfFlag
    %fullname=sprintf('%s/%s.%s',printoutDir,name,'eps');
    %fulln=sprintf('%s/%s',printoutDir,name);
    if verboseflag
        fprintf('printing to %s.pdf \n',fulln);
    end
    export_fig(fulln,'-pdf','-transparent')
end


%set(gcf, 'PaperPositionMode', 'auto');
%print(gcf,'-depsc2',fullname);
%fixPSlinestyle(fullname);
%exportfig(figHandle,fullname,'color','rgb','width',20);

%saveas(figHandle,fulln,'pdf')

%export_fig(fullname,'-png')
%savefigure(fullname,'fontsize',14,'format','eps');
%fixPSlinestyle(fullname);
%epswrite(fullname,'size','screen','units','centimeters','boundingbox','tight');


if figFlag
    
    ffDir=[printoutDir '/figFiles'];
    dirStat=exist(ffDir,'dir');
    if dirStat==0
        mkdir(ffDir)
        dirStat=exist(ffDir,'dir');
        
        if dirStat~=7
            error('PRINTOUT_FIG - something wrong with figFiles directory: %s',printoutDir);
        end
    end
    fullname=sprintf('%s/%s.%s',ffDir,name,'fig');
    saveas(figHandle,fullname,'fig')
end

