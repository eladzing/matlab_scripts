function galMovie_plot_gasMass(zr,comp,indx,varargin )
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
semilogx(zr,comp.sfg,greyType,'color',grey,'linewidth',1);
hold on
semilogx(zr,comp.cdn,greyType,'color',grey,'linewidth',1);
semilogx(zr,comp.cdi,greyType,'color',grey,'linewidth',1);
semilogx(zr,comp.wrm,greyType,'color',grey,'linewidth',1);
semilogx(zr,comp.hot,greyType,'color',grey,'linewidth',1);

%plot history in color
h(1)=semilogx(zr(indx:end),comp.sfg(indx:end),'-','color',colorSet(4,:),'linewidth',1.5,'DisplayName','SF');
h(2)=semilogx(zr(indx:end),comp.cdn(indx:end),'-','color',colorSet(2,:),'linewidth',1.5,'DisplayName','Cold Dense');
h(3)=semilogx(zr(indx:end),comp.cdi(indx:end),'-','color',colorSet(3,:),'linewidth',1.5,'DisplayName','Cold Dilute');
h(4)=semilogx(zr(indx:end),comp.wrm(indx:end),'-','color',colorSet(5,:),'linewidth',1.5,'DisplayName','Warm');
h(5)=semilogx(zr(indx:end),comp.hot(indx:end),'-','color',colorSet(1,:),'linewidth',1.5,'DisplayName','Hot');

% plot
semilogx(zr(indx),comp.sfg(indx),'d','MarkerFaceColor',colorSet(4,:),'markersize',8)
semilogx(zr(indx),comp.cdn(indx),'d','MarkerFaceColor',colorSet(2,:),'markersize',8)
semilogx(zr(indx),comp.cdi(indx),'d','MarkerFaceColor',colorSet(3,:),'markersize',8)
semilogx(zr(indx),comp.wrm(indx),'d','MarkerFaceColor',colorSet(5,:),'markersize',8)
semilogx(zr(indx),comp.hot(indx),'d','MarkerFaceColor',colorSet(1,:),'markersize',8)

if ylabelLeftFlag
ylabelmine('Mass Fraction')
end

yyaxis right

h(6)=semilogx(zr,comp.gmass./1e10,'k-.','linewidth',0.8,'DisplayName','gas mass');%,'color',grey)
hold on
semilogx(zr(indx),comp.gmass(indx)./1e10,'dk',...
    'markersize',10);%,'MarkerFaceColor',colorSet(3,:));

if ylabelRightFlag
ylabelmine('Gas Mass $[\mathrm{10^{10}M_\odot}]$')
end

if legendFlag
    hl=legend(h);
    set(hl,'Interpreter','latex','Fontsize',10);
end

if xlabelFlag
xlabelmine('$\log(1+z)$')
end

%set(gca,'Fontsize',12)



end

