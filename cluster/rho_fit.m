function fits= rho_fit(lrp,lrot,lrodm,lrog,rv)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(LRP,LROT,LRODM,LROG)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  3
%   Number of fits:  8


% Data from dataset "lrot vs. lrp":
%    X = lrp:
%    Y = lrot:
%    Unweighted

% Data from dataset "lrodm vs. lrp":
%    X = lrp:
%    Y = lrodm:
%    Unweighted

% Data from dataset "lrog vs. lrp":
%    X = lrp:
%    Y = lrog:
%    Unweighted
%
% This function was automatically generated on 16-Mar-2010 16:22:29


% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[577 56 678 566]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
grid(ax_,'on');
axes(ax_); hold on;


% --- Plot data originally in dataset "lrot vs. lrp"
lrp = lrp(:);
lrot = lrot(:);
h_ = line(lrp,lrot,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(lrp));
xlim_(2) = max(xlim_(2),max(lrp));
legh_(end+1) = h_;
legt_{end+1} = 'lrot vs. lrp';

% --- Plot data originally in dataset "lrodm vs. lrp"
lrodm = lrodm(:);
h_ = line(lrp,lrodm,'Parent',ax_,'Color',[0.333333 0.666667 0],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(lrp));
xlim_(2) = max(xlim_(2),max(lrp));
legh_(end+1) = h_;
legt_{end+1} = 'lrodm vs. lrp';

% --- Plot data originally in dataset "lrog vs. lrp"
lrog = lrog(:);
h_ = line(lrp,lrog,'Parent',ax_,'Color',[0 0 0],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(lrp));
xlim_(2) = max(xlim_(2),max(lrp));
legh_(end+1) = h_;
legt_{end+1} = 'lrog vs. lrp';


% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-2.2834115786954694016, 0.78409290497960859589]);
end

iif=0;

% --- Create fit "fit 1"

% Apply exclusion rule "10cellmin"
rmax=log10(rv');  %max([0.5 log10(rv)]);
%rmax=0.5;
ex_ = (lrp <= -1.25 | lrp >= rmax);
ok_ = isfinite(lrp) & isfinite(lrot);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [15 0.10000000000000000555 ];
ft_ = fittype('ros-(x-log10(rs)+2*log10(1+10^(x-log10(rs))))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'ros', 'rs'});

% Fit this model using new data
if sum(~ex_(ok_))<2  %% too many points excluded
    error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 1','10cellmin')
else
    [cf_ gof] = fit(lrp(ok_),lrot(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
end

% Or use coefficients from the original fit:
if 0
    cv_ = { 15.330523474061873657, 0.38621362822991411878};
    [cf_ gof] = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'nfw';

iif=iif+1;
fits(iif).data='tot';
fits(iif).model='nfw';
fits(iif).cft=cf_;
fits(iif).gof=gof;

% % --- Create fit "fit 2"
% 
% % Apply exclusion rule "10cellmin"
% ex_ = (lrp <= -1.25 | lrp >= rmax);
% ok_ = isfinite(lrp) & isfinite(lrot);
% if ~all( ok_ )
%     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
%         'Ignoring NaNs and Infs in data' );
% end
% st_ = [1 15 0.10000000000000000555 ];
% ft_ = fittype('ros-(a*(x-log10(rs))+(3-a)*log10(1+10^(x-log10(rs))))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a', 'ros', 'rs'});
% 
% % Fit this model using new data
% if sum(~ex_(ok_))<2  %% too many points excluded
%     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 2','10cellmin')
% else
%     [cf_ gof] = fit(lrp(ok_),lrot(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% end
% 
% % Or use coefficients from the original fit:
% if 0
%     cv_ = { 1.4724947447409753032, 14.484590198014082318, 0.79891685214253593816};
%     [cf_ gof] = cfit(ft_,cv_{:});
% end
% 
% % Plot this fit
% h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[0 0 1],...
%     'LineStyle','-', 'LineWidth',2,...
%     'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'nfw-a';
% 
% iif=iif+1;
% fits(iif).data='tot';
% fits(iif).model='nfw-a';
% fits(iif).cft=cf_;
% fits(iif).gof=gof;
% 
% % % --- Create fit "fit 3"
% % 
% % % Apply exclusion rule "10cellmin"
% % ex_ = (lrp <= -1.25 | lrp >= rmax);
% % ok_ = isfinite(lrp) & isfinite(lrot);
% % if ~all( ok_ )
% %     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
% %         'Ignoring NaNs and Infs in data' );
% % end
% % st_ = [2 15 0.10000000000000000555 ];
% % ft_ = fittype('ros-a*(x-log10(rs))',...
% %     'dependent',{'y'},'independent',{'x'},...
% %     'coefficients',{'a', 'ros', 'rs'});
% % 
% % % Fit this model using new data
% % if sum(~ex_(ok_))<2  %% too many points excluded
% %     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 3','10cellmin')
% % else
% %     [cf_ gof] = fit(lrp(ok_),lrot(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% % end
% % 
% % % Or use coefficients from the original fit:
% % if 0
% %     cv_ = { 2.0323704220951901434, 15.269302868699702103, 0.17991535674742747952};
% %     [cf_ gof] = cfit(ft_,cv_{:});
% % end
% % 
% % % Plot this fit
% % h_ = plot(cf_,'fit',0.95);
% % legend off;  % turn off legend from plot method call
% % set(h_(1),'Color',[0.666667 0.333333 0],...
% %     'LineStyle','-', 'LineWidth',2,...
% %     'Marker','none', 'MarkerSize',6);
% % legh_(end+1) = h_(1);
% % legt_{end+1} = 'p-law';
% % 
% % iif=iif+1;
% % fits(iif).data='tot';
% % fits(iif).model='p-law';
% % fits(iif).cft=cf_;
% % fits(iif).gof=gof;
% 
% % --- Create fit "fit 4"
% 
% % Apply exclusion rule "10cellmin"
% ex_ = (lrp <= -1.25 | lrp >= rmax);
% ok_ = isfinite(lrp) & isfinite(lrodm);
% if ~all( ok_ )
%     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
%         'Ignoring NaNs and Infs in data' );
% end
% st_ = [15 0.10000000000000000555 ];
% ft_ = fittype('ros-(x-log10(rs)+2*log10(1+10^(x-log10(rs))))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'ros', 'rs'});
% 
% % Fit this model using new data
% if sum(~ex_(ok_))<2  %% too many points excluded
%     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 4','10cellmin')
% else
%     [cf_ gof] = fit(lrp(ok_),lrodm(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% end
% 
% % Or use coefficients from the original fit:
% if 0
%     cv_ = { 14.841878077007091719, 0.36681638276236383511};
%     [cf_ gof] = cfit(ft_,cv_{:});
% end
% 
% 
% % Plot this fit
% h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[0.333333 0.333333 0.333333],...
%     'LineStyle','-', 'LineWidth',2,...
%     'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'nfw';
% 
% iif=iif+1;
% fits(iif).data='dm';
% fits(iif).model='nfw';
% fits(iif).cft=cf_;
% fits(iif).gof=gof;
% 
% % --- Create fit "fit 5"
% 
% % Apply exclusion rule "10cellmin"
% ex_ = (lrp <= -1.25 | lrp >= rmax);
% ok_ = isfinite(lrp) & isfinite(lrodm);
% if ~all( ok_ )
%     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
%         'Ignoring NaNs and Infs in data' );
% end
% st_ = [1 15 0.10000000000000000555 ];
% ft_ = fittype('ros-(a*(x-log10(rs))+(3-a)*log10(1+10^(x-log10(rs))))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a', 'ros', 'rs'});
% 
% % Fit this model using new data
% if sum(~ex_(ok_))<2  %% too many points excluded
%     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 5','10cellmin')
% else
%     [cf_ gof] = fit(lrp(ok_),lrodm(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% end
% 
% % Or use coefficients from the original fit:
% if 0
%     cv_ = { 1.5264466642447322986, 13.880186690822363005, 0.83672889441757369866};
%     [cf_ gof] = cfit(ft_,cv_{:});
% end
% 
% % Plot this fit
% h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[1 0 1],...
%     'LineStyle','-', 'LineWidth',2,...
%     'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'nfw-a';
% 
% 
% iif=iif+1;
% fits(iif).data='dm';
% fits(iif).model='nfw-a';
% fits(iif).cft=cf_;
% fits(iif).gof=gof;
% 
% % % --- Create fit "fit 6"
% % 
% % % Apply exclusion rule "10cellmin"
% % ex_ = (lrp <= -1.25 | lrp >= rmax);
% % ok_ = isfinite(lrp) & isfinite(lrodm);
% % if ~all( ok_ )
% %     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
% %         'Ignoring NaNs and Infs in data' );
% % end
% % st_ = [2 15 0.10000000000000000555 ];
% % ft_ = fittype('ros-a*(x-log10(rs))',...
% %     'dependent',{'y'},'independent',{'x'},...
% %     'coefficients',{'a', 'ros', 'rs'});
% % 
% % % Fit this model using new data
% % if sum(~ex_(ok_))<2  %% too many points excluded
% %     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 6','10cellmin')
% % else
% %     [cf_ gof] = fit(lrp(ok_),lrodm(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% % end
% % 
% % % Or use coefficients from the original fit:
% % if 0
% %     cv_ = { 2.0533536145838984233, 15.018550609167032661, 0.13204243206550630996};
% %     [cf_ gof] = cfit(ft_,cv_{:});
% % end
% % 
% % h_ = plot(cf_,'fit',0.95);
% % legend off;  % turn off legend from plot method call
% % set(h_(1),'Color',[1 1 0],...
% %     'LineStyle','-', 'LineWidth',2,...
% %     'Marker','none', 'MarkerSize',6);
% % legh_(end+1) = h_(1);
% % legt_{end+1} = 'p-law';
% % 
% % 
% % iif=iif+1;
% % fits(iif).data='dm';
% % fits(iif).model='p-law';
% % fits(iif).cft=cf_;
% % fits(iif).gof=gof;
% 
% % --- Create fit "fit 7"
% 
% % Apply exclusion rule "10cellmin"
% ex_ = (lrp <= -1.25 | lrp >= rmax);
% ok_ = isfinite(lrp) & isfinite(lrog);
% if ~all( ok_ )
%     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
%         'Ignoring NaNs and Infs in data' );
% end
% st_ = [15 0.10000000000000000555 ];
% ft_ = fittype('ros-log10(1+10^(2*(x-log10(rs))))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'ros', 'rs'});
% 
% % Fit this model using new data
% if sum(~ex_(ok_))<2  %% too many points excluded
%     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 7','10cellmin')
% else
%     [cf_ gof] = fit(lrp(ok_),lrog(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% end
% 
% % Or use coefficients from the original fit:
% if 0
%     cv_ = { 14.59522437733601663, 0.11928766388890452155};
%     [cf_ gof] = cfit(ft_,cv_{:});
% end
% 
% % Plot this fit
% h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[1 0.666667 0.333333],...
%     'LineStyle','-', 'LineWidth',2,...
%     'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'isothermal-core';
% 
% 
% iif=iif+1;
% fits(iif).data='gas';
% fits(iif).model='iso-core';
% fits(iif).cft=cf_;
% fits(iif).gof=gof;
% 
% % --- Create fit "fit 8"
% 
% % Apply exclusion rule "10cellmin"
% ex_ = (lrp <= -1.25 | lrp >= rmax);
% ok_ = isfinite(lrp) & isfinite(lrog);
% if ~all( ok_ )
%     warning( 'GenerateMFile:IgnoringNansAndInfs', ...
%         'Ignoring NaNs and Infs in data' );
% end
% st_ = [2 15 0.10000000000000000555 ];
% ft_ = fittype('ros-log10(1+10^(a*(x-log10(rs))))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a', 'ros', 'rs'});
% 
% % Fit this model using new data
% if sum(~ex_(ok_))<2  %% too many points excluded
%     error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 8','10cellmin')
% else
%     [cf_ gof] = fit(lrp(ok_),lrog(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
% end
% 
% % Or use coefficients from the original fit:
% if 0
%     cv_ = { 2.2159999767393379067, 14.478982580926221146, 0.16455715720390756696};
%     [cf_ gof] = cfit(ft_,cv_{:});
% end
% 
% h_ = plot(cf_,'fit',0.95);
% legend off;  % turn off legend from plot method call
% set(h_(1),'Color',[0.666667 0.666667 0.666667],...
%     'LineStyle','-', 'LineWidth',2,...
%     'Marker','none', 'MarkerSize',6);
% legh_(end+1) = h_(1);
% legt_{end+1} = 'p-law-core';
% 
% 
% iif=iif+1;
% fits(iif).data='gas';
% fits(iif).model='p-law-core';
% fits(iif).cft=cf_;
% fits(iif).gof=gof;
% 
% yl=[min(lrog) 0.5*(min(lrog)+max(lrot))];
% lrv=[log10(rv) log10(rv)];
% h_=line(lrv,yl,'Color','k','LineStyle','--');
% legh_(end+1) = h_;
% legt_{end+1} = 'R_{vir}';
% 
% % Done plotting data and fits.  Now finish up loose ends.
% hold off;
% leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
% h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
% set(h_,'Interpreter','none');
% xlabel(ax_,'log_{10}(r) [Mpc]');               % remove x label
% ylabel(ax_,'\rho [M_{sun}/Mpc^3]');               % remove y label
