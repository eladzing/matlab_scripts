function plot_param_vs_radius_2Dhistogram(limRad,limParam,bird,varargin )
%PLOT_PARAM_VS_2DHISTOGRAM plot a 2D histogram (bird plot) where one of the
%paramters is the radial distance from the center of the object
%

% defualts
bartag='Mass';
bmap='*Spectral';
logFlag=true;
xlab="radial distance ";
ylab="";
%histFlag=false;
titletag='';
%radNorm=1;
%guideFlag=false;
%tmpRoFlag=false;
colLim=[1 1];
%hf=figure('Visible','off');
axesHandle=0;

plotBird=bird'; %% default - show what came in 


% guide(1).tmp=[7 8];
% guide(1).ro=[0 NaN];
% guide(1).a=-1;
% guide(1).c=guide(1).ro(1)-guide(1).a.*guide(1).tmp(1);
% guide(1).ro(2)=guide(1).a.*guide(1).tmp(2)+guide(1).c;
% 
% guide(2).tmp=[4 5];
% guide(2).ro=[-5 NaN];
% guide(2).a=1.5;
% guide(2).c=guide(2).ro(1)-guide(2).a.*guide(2).tmp(1);
% guide(2).ro(2)=guide(2).a.*guide(2).tmp(2)+guide(2).c;


i=1;

while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'brewermap','brewer'}
            i=i+1;
            bmap=varargin{i};
        case {'cmap','colormap'}
            i=i+1;
            cmap=varargin{i};
        case {'tag','bartag','bartitle','bar'}
            i=i+1;
            bartag=varargin{i};
        case {'title','tit'}
            i=i+1;
            titletag=varargin{i};
        case {'xlabel','xlab','rlabel','rlab'}
            i=i+1;
            xlab=varargin{i};
        case {'ylabel','ylab','parlabel','parlab','paramlabel','paramlab'}
            i=i+1;
            ylab=varargin{i};
        case 'nolog'
            logFlag=false;
        case 'log'
            logFlag=true;
        % case {'guide','guideline','lines'}
        %     guideFlag=true;
        case {'caxis','clims','clim'}
            i=i+1;
            colLim=varargin{i};
        % case{'tro','t-ro','tn','t-n'}
        %     tmpRoFlag=true;
        % case{'rot','ro-t','nt','n-t'}
        %     tmpRoFlag=false;
        case{'handle','fighandle','fig','figure'}
            i=i+1;
            hf=varargin{i};
        case{'axes','axeshandle'}
            i=i+1;
            axesHandle=varargin{i};
        case{'hist','histogram','showhist'}
            histFlag=true;
        case{'fractional','frac','percent','perc'}
            plotBird=bird./sum(sum(bird)).*100;
            
               bartag='$ \%$ Mass';
       


        otherwise
    end
    i=i+1;
end

if logFlag
    plotBird=log10(plotBird);
end

if ~exist('cmap','var')
    cmap=brewermap(256,bmap);   
end
cmap(1,:)=[1 1 1];
if ~exist('hf')
    hf=myFigure;
    %hf.visibility='on';
else
    figure(hf);
    %visState=hf1.Visible;
    %figure(hf);%,'Visible',visState)
    %hf=hf1;
    %hf.Visible=visState;
end

if axesHandle~=0
    axes(axesHandle);
end

if isempty(ylab)
    warning('%s - no y-label given, How will people know what you are plotting?',current_function().upper);
    ylab="unknown parameter";
end

imsc(limRad,limParam,plotBird,cmap);

set(gca,'Ydir','normal','Fontsize',12)

if colLim==[1 1]
    colLim(1)=min(min(plotBird));
    
    colLim(2)=max(max(plotBird));
end


% if guideFlag
%     hold on
%     for i=1:2
%         if tmpRoFlag
%             plot(guide(i).tmp,guide(i).ro,'--k')
%         else
%             plot(guide(i).ro,guide(i).tmp,'--k')
%         end
%     end
% end



grid('on')

colormap(cmap);
hb=colorbar;

barTitle(hb,bartag)
clim(colLim)

if string(xlab)~="none"
    xlabelmine(xlab);
end
if string(ylab)~="none"
    ylabelmine(ylab);
end

if ~isempty(titletag)
    titlemine(titletag);
end
myAxis;

hf1=hf;
%% add single value histogram on the side
% if histFlag  % show histogram
%
%     histT=sum(bird,1);
%     histRo=sum(bird,2);
%
%     sz=size(bird)-1;
%     axisT=limT(1):diff(limT)/sz(1):limT(2);
%     axisRo=limRo(1):diff(limRo)/sz(2):limRo(2);
%
%     
%     
%     bar(axisRo,log10(histRo),'b')
%     xlim(limRo)
%     barh(axisT,log10(histT),'k')
%     ylim(limT)
    
%end 




end

