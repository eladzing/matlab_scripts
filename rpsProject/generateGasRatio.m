function res = generateGasRatio(stellarMass,varargin)
% GENERATEGASRATIO - Generate the gas to stellar mass ratio for a galaxy
%   Using the distributions from TNG100/TNG300, generate random values of
%   the gas to stellar mass ratio.
%   stellarMass - stellar mass of the galaxy
%   TNG100/TNG300 - Use distribution from simulations

simType='TNG100';
snap=99;
histBin=200;
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        case{'100',100,'tng100'}
            simType='TNG100';
        case{'300',300,'tng300'}
            simType='TNG300';
        case{'snap'}
            i=i+1;
            snap=varargin{i};
            if ~any(snap==[50,67,99])
                error('GENERATEGASRATIO - Fit parameteres do not exist for snapshot: %i',num2str(snap));
            end
        case{'histbin'}
            i=i+1;
            histBin=varargin{i};
        otherwise
            error('GENERATEGASRATIO - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end

%% read in gaussian fit files

ft=read_fgs_gaussian_fit_files(simType,snap,histBin);



%% run over mass bins
binInd=discretize(log10(stellarMass),ft.binEdge);

resTemp=zeros(size(stellarMass));

for i=1:length(ft.binEdge)-1
    mask=binInd==i;
    
    
    if any(mask)
        %% create pdf & cdf
        fgs=-5:0.001:0.5;
        
        pdf=ft.a1(i).*exp(-0.5.*((fgs-ft.mu1(i))./ft.sig1(i)).^2) + ...
            ft.a2(i).*exp(-0.5.*((fgs-ft.mu2(i))./ft.sig2(i)).^2);
        
        
        cdf=cumtrapz(fgs,pdf)./trapz(fgs,pdf);
        
        
        %     figure
        %     yyaxis left
        %     plot(fgs,pdf)
        %     ylabelmine('PDF');
        %     yyaxis right
        %     plot(fgs,cdf)
        %     ylabelmine('CDF');
        
        
        
        %% generate random values
        [cdf, cdfMask] = unique(cdf);
        fgs=fgs(cdfMask);
        
        rv=rand(1,sum(mask));
        %fgsGen=interp1(cdf,fgs,rv);
        
        resTemp(mask)=interp1(cdf,fgs,rv);
        
    end
    
    
end
res=10.^resTemp;
%res=resTemp;

% resTemp=resTemp(resTemp>0);
%
% df=nOrb-length(resTemp);
% cnt=0;
% while df>0 || cnt>20
%     cnt=cnt+1;
%
%     rv=rand(1,df);
%     resT=interp1(cdf,vv,rv);
%     resT=resT(resT>0);
%
%     resTemp=cat(2,resTemp,resT);
%
%     df=nOrb-length(resTemp);
% end
%






