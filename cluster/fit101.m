function createFit(rp,rot,lrp,lrot,rodm,lrodm)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(RP,ROT,LRP,LROT,RODM,LRODM)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  4
%   Number of fits:  4


% Data from dataset "rot vs. rp":
%    X = rp:
%    Y = rot:
%    Unweighted

% Data from dataset "lrot vs. lrp":
%    X = lrp:
%    Y = lrot:
%    Unweighted

% Data from dataset "rodm vs. rp":
%    X = rp:
%    Y = rodm:
%    Unweighted

% Data from dataset "lrodm vs. lrp":
%    X = lrp:
%    Y = lrodm:
%    Unweighted
%
% This function was automatically generated on 02-Mar-2010 16:39:24

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[552 144 725 483]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 .5 1 .5]);
ax2_ = axes;
set(ax2_,'Units','normalized','OuterPosition',[0 0 1 .5]);
set(ax2_,'Box','on');
legrh_ = []; legrt_ = {};
set(ax_,'Box','on');
axes(ax_); hold on;


% --- Plot data originally in dataset "rot vs. rp"
rp = rp(:);
rot = rot(:);
h_ = line(rp,rot,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(rp));
xlim_(2) = max(xlim_(2),max(rp));
legh_(end+1) = h_;
legt_{end+1} = 'rot vs. rp';

% --- Plot data originally in dataset "lrot vs. lrp"
lrp = lrp(:);
lrot = lrot(:);
h_ = line(lrp,lrot,'Parent',ax_,'Color',[0.333333 0.666667 0],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(lrp));
xlim_(2) = max(xlim_(2),max(lrp));
legh_(end+1) = h_;
legt_{end+1} = 'lrot vs. lrp';

% --- Plot data originally in dataset "rodm vs. rp"
rodm = rodm(:);
% This dataset does not appear on the plot
% Add it to the plot by removing the if/end statements that follow
% and by selecting the desired color and marker
if 0
    h_ = line(rp,rodm,'Color','r','Marker','.','LineStyle','none');
    xlim_(1) = min(xlim_(1),min(rp));
    xlim_(2) = max(xlim_(2),max(rp));
    legh_(end+1) = h_;
    legt_{end+1} = 'rodm vs. rp';
end       % end of "if 0"

% --- Plot data originally in dataset "lrodm vs. lrp"
lrodm = lrodm(:);
% This dataset does not appear on the plot
% Add it to the plot by removing the if/end statements that follow
% and by selecting the desired color and marker
if 0
    h_ = line(lrp,lrodm,'Color','r','Marker','.','LineStyle','none');
    xlim_(1) = min(xlim_(1),min(lrp));
    xlim_(2) = max(xlim_(2),max(lrp));
    legh_(end+1) = h_;
    legt_{end+1} = 'lrodm vs. lrp';
end       % end of "if 0"

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
    set(ax2_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-2.283436761663227621, 0.78663638472317520822]);
    set(ax2_,'XLim',[-2.283436761663227621, 0.78663638472317520822]);
end


% --- Create fit "fit 1"

% Apply exclusion rule "2"
ex_ = (rp <= -1.5 | rp >= 0.5);
ok_ = isfinite(rp) & isfinite(rot);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [15 0.10000000000000000555 ];
ft_ = fittype('ros-(x-log10(rs))-2*log10(1+10^(x-log10(rs)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'ros', 'rs'});

% Fit this model using new data
if sum(~ex_(ok_))<2  %% too many points excluded
    error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 1','2')
else
    cf_ = fit(rp(ok_),rot(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
end

% Or use coefficients from the original fit:
if 0
    cv_ = { 15.184640155120312954, 0.44984073871217195029};
    cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';
res_ = rot(~ex_) - cf_(rp(~ex_));
[x_,i_] = sort(rp(~ex_));
axes(ax2_); hold on;
h_ = line(x_,res_(i_),'Parent',ax2_,'Color',[1 0 0],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',6);
axes(ax_); hold on;
legrh_(end+1) = h_;
legrt_{end+1} = 'fit 1';

% --- Create fit "fit 2"

% Apply exclusion rule "2"
ex_ = (lrp <= -1.5 | lrp >= 0.5);
ok_ = isfinite(lrp) & isfinite(lrot);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [15 0.10000000000000000555 ];
ft_ = fittype('ros-(x-log10(rs))-2*log10(1+10^(x-log10(rs)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'ros', 'rs'});

% Fit this model using new data
if sum(~ex_(ok_))<2  %% too many points excluded
    error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 2','2')
else
    cf_ = fit(lrp(ok_),lrot(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
end

% Or use coefficients from the original fit:
if 0
    cv_ = { 15.388617599584742734, 0.36445549966839579925};
    cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[0 0 1],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 2';
res_ = lrot(~ex_) - cf_(lrp(~ex_));
[x_,i_] = sort(lrp(~ex_));
axes(ax2_); hold on;
h_ = line(x_,res_(i_),'Parent',ax2_,'Color',[0 0 1],...
    'LineStyle','-', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',6);
axes(ax_); hold on;
legrh_(end+1) = h_;
legrt_{end+1} = 'fit 2';

% --- Create fit "fit 3"

% Apply exclusion rule "2"
ex_ = (rp <= -1.5 | rp >= 0.5);
ok_ = isfinite(rp) & isfinite(rodm);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [15 0.10000000000000000555 ];
ft_ = fittype('ros-(x-log10(rs))-2*log10(1+10^(x-log10(rs)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'ros', 'rs'});

% Fit this model using new data
if sum(~ex_(ok_))<2  %% too many points excluded
    error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 3','2')
else
    cf_ = fit(rp(ok_),rodm(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
end

% Or use coefficients from the original fit:
if 0
    cv_ = { 14.183554501964650285, 0.4325069333346459044};
    cf_ = cfit(ft_,cv_{:});
end

% This fit does not appear on the plot

% --- Create fit "fit 4"

% Apply exclusion rule "2"
ex_ = (lrp <= -1.5 | lrp >= 0.5);
ok_ = isfinite(lrp) & isfinite(lrodm);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [15 0.10000000000000000555 ];
ft_ = fittype('ros-(x-log10(rs))-2*log10(1+10^(x-log10(rs)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'ros', 'rs'});

% Fit this model using new data
if sum(~ex_(ok_))<2  %% too many points excluded
    error('Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.','fit 4','2')
else
    cf_ = fit(lrp(ok_),lrodm(ok_),ft_,'Startpoint',st_,'Exclude',ex_(ok_));
end

% Or use coefficients from the original fit:
if 0
    cv_ = { 14.40720597022787075, 0.3438177591146819867};
    cf_ = cfit(ft_,cv_{:});
end

% This fit does not appear on the plot

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
h_ = legend(ax2_,legrh_,legrt_,leginfo_{:});
set(h_,'Interpreter','none');
xlabel(ax2_,'');
ylabel(ax2_,'');
title(ax_,'Data and Fits');
title(ax2_,'Residuals');
