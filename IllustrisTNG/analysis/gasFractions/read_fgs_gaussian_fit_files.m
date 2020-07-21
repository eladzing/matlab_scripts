function res = read_fgs_gaussian_fit_files(simType,snap,histBin)
% READ_FGS_GAUSSIAN_FIT_FILES - read in the fit files and return the fit
% parameters
%   Detailed explanation goes here

%simType='TNG100';
%snap=99;
%histBin=200;

global DEFAULT_MATFILE_DIR
if contains(DEFAULT_MATFILE_DIR,'IllustrisTNG') 
    matDir=DEFAULT_MATFILE_DIR;
elseif contains(lower(DEFAULT_MATFILE_DIR),'rps')
    global DEFAULT_TNG_MATFILE_DIR
    matDir=DEFAULT_TNG_MATFILE_DIR;
else 
    error('READ_FGS_GAUSSIAN_FIT_FILES - incorrect MATFILE_DIR: %s',DEFAULT_MATFILE_DIR)
end
    


load([matDir '/fgs_stellarMass_gaussianFit_byBin_snp_' num2str(snap) '_TNG100.mat'],'ft')
ft1=ft;

load([matDir '/fgs_stellarMass_gaussianFit_byBin_snp_' num2str(snap) '_TNG300.mat'],'ft')
ft3=ft;

fitInd=floor(histBin/25/2+1);
if ft1(fitInd).histBinNum~=histBin || ft3(fitInd).histBinNum~=histBin
    error('GENERATEGASRATIO - Something wrong with finding correct fit index');
end
binEdge=ft3(fitInd).binEdges;

 % We use TNG300 distributions for all stellar masses above 1e10.75; 

 mu1= ft3(fitInd).mu1;
 sig1= ft3(fitInd).sig1;
 mu2= ft3(fitInd).mu2;
 sig2= ft3(fitInd).sig2;
 a1= ft3(fitInd).a1;
 a2= ft3(fitInd).a2;

if strcmp(simType,'TNG100')
        bins=binEdge<=10.75;
        % We use TNG300 distributions for all stellar masses above 1e11; 
        
        mu1(bins)= ft1(fitInd).mu1(bins);
        sig1(bins)= ft1(fitInd).sig1(bins);
        mu2(bins)= ft1(fitInd).mu2(bins);
        sig2(bins)= ft1(fitInd).sig2(bins);
        a1(bins)= ft1(fitInd).a1(bins);
        a2(bins)= ft1(fitInd).a2(bins);
end



res.mu1=mu1;
res.sig1= sig1;
res.mu2=mu2;
res.sig2=sig2;
res.a1=a1;
res.a2=a2;
res.binEdge=binEdge;



end

