function galMovie_plot_gasEnergy(zr,zr2,comp,ii,varargin )
%GALMOVIE_PLOT_GASMASS auxillary function to make a subplot in the movie galaxy analyssi movie
%   Detailed explanation goes here

%% perliminaries
colorSet=brewermap(9,'Set1');
grey=colorSet(9,:);
greyType=':';

legendFlag=false;
xlabelFlag=false;
ylabelLeftFlag=false;
ylabelRightFlag=false;

%% parse flags 
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'legend'}
            legendFlag=true;
        case{'xlabel','xlab','xl'}
        xlabelFlag=true;
        
        case{'ylabel','ylab','ylabelleft','yleft','yl'}
        ylabelLeftFlag=true;
        
        case{'ylabelright','ylabright','yright','yr'}
        ylabelRightFlag=true;
        otherwise
            error(' galMovie_plot_gasMass - Illegal Argument: %s',varargin{i});
    end
    i=i+1;
end


%% plot
h=[];

yyaxis left

% plot grey background
loglog(zr,comp.sfg,greyType,'color',grey,'linewidth',1);
hold on
loglog(zr,comp.cdn,greyType,'color',grey,'linewidth',1);
loglog(zr,comp.cdi,greyType,'color',grey,'linewidth',1);
loglog(zr,comp.wrm,greyType,'color',grey,'linewidth',1);
loglog(zr,comp.hot,greyType,'color',grey,'linewidth',1);

%plot history in color
h(1)=loglog(zr(ii:end),comp.sfg(ii:end),'-','color',colorSet(4,:),'linewidth',1.5,'DisplayName','SF');
h(2)=loglog(zr(ii:end),comp.cdn(ii:end),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','Cold Dense');
h(3)=loglog(zr(ii:end),comp.cdi(ii:end),'-','color',colorSet(3,:),'linewidth',1.5,'DisplayName','Cold Dilute');
h(4)=loglog(zr(ii:end),comp.wrm(ii:end),'-','color',colorSet(5,:),'linewidth',1.5,'DisplayName','Warm');
h(5)=loglog(zr(ii:end),comp.hot(ii:end),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','Hot');

indx0=find(zr2==zr(ii));
if ~isempty(indx0)
    indx=indx0;
end


% plot
loglog(zr2(indx:end),comp.sfg2(indx:end),'d','MarkerFaceColor',colorSet(4,:),'markersize',8,'color',colorSet(4,:))
loglog(zr2(indx:end),comp.cdn2(indx:end),'d','MarkerFaceColor',colorSet(2,:),'markersize',8,'color',colorSet(2,:))
loglog(zr2(indx:end),comp.cdi2(indx:end),'d','MarkerFaceColor',colorSet(3,:),'markersize',8,'color',colorSet(3,:))
loglog(zr2(indx:end),comp.wrm2(indx:end),'d','MarkerFaceColor',colorSet(5,:),'markersize',8,'color',colorSet(5,:))
loglog(zr2(indx:end),comp.hot2(indx:end),'d','MarkerFaceColor',colorSet(1,:),'markersize',8,'color',colorSet(1,:))


grid 
%grid minor


if ylabelLeftFlag
ylabelmine('Mass Fraction')
end

% yyaxis right
% % 
% h(6)=loglog(zr,comp.gmass./1e10,'k-.','linewidth',0.8,'DisplayName','gas mass');%,'color',grey)
% hold on
% loglog(zr(indx),comp.gmass(indx)./1e10,'dk',...
%     'markersize',10);%,'MarkerFaceColor',colorSet(3,:));

if ylabelRightFlag
%ylabelmine('Gas Mass $[\mathrm{10^{10}M_\odot}]$')
ylabelmine('$E_{\mathrm{dis}}/ M_{\mathrm{g}}\,[\mathrm{erg\,yr^{-1}\,M_\odot}]$')

end

if legendFlag
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',10);
end

if xlabelFlag
xlabelmine('$\log(1+z)$')
end

set(gca,'Fontsize',12)



end

