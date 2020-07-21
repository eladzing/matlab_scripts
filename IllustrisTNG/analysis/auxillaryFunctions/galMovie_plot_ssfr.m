function  galMovie_plot_ssfr(zr,ghs,indx,varargin )  
%GALMOVIE_PLOT_ssfr auxillary function to make a subplot in the movie galaxy analyssi movie
%   Detailed explanation goes here


%%     obsolete 
%%


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
semilogy(zr,ghs.ssfr,greyType,'color',grey)
hold on
semilogy(zr(indx:end),ghs.ssfr(indx:end),'-','color',colorSet(2,:),'linewidth',1.5)
semilogy(zr(indx),ghs.ssfr(indx),'d','color',colorSet(2,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(2,:));
if ylabelLeftFlag
ylabelmine('ssfr $[\mathrm{yr^{-1}}]$')
end

ylim([1e-17 1e-8])


yyaxis right

semilogy(zr,ghs.stellarMass./1e10,greyType,'color',grey) 
hold on
h(1)=semilogy(zr(indx:end),ghs.stellarMass(indx:end)./1e10,'-',...
    'color',colorSet(1,:),'linewidth',1.5,'DisplayName','Stellar');
semilogy(zr(indx),ghs.stellarMass(indx)./1e10,'d','color',colorSet(1,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(1,:));
    

semilogy(zr,ghs.inGal.gasMass(2,:)./1e10,greyType,'color',grey) 
hold on
h(2)=semilogy(zr(indx:end),ghs.inGal.gasMass(2,indx:end)./1e10,'-',...
    'color',colorSet(3,:),'linewidth',1.5,'DisplayName','Gas');
semilogy(zr(indx),ghs.inGal.gasMass(2,indx)./1e10,'d','color',colorSet(3,:),...
    'linewidth',1.5,'markersize',10,'MarkerFaceColor',colorSet(3,:));



if ylabelRightFlag
ylabelmine('Gas / Stellar Mass $[\mathrm{10^{10}M_\odot}]$')
end


if legendFlag
hl=legend(h);
set(hl,'Interpreter','latex','Fontsize',14);
end


if xlabelFlag
    xlabelmine('redshift')
end


set(gca,'Fontsize',12)



end

