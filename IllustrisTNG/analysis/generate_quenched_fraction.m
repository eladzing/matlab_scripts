function res =generate_quenched_fraction(fofs,subs,varargin)
%GENERATE_QUENCHED_FRACTION prepare quenched fraction profiles  binned to
%host and stellar mass
%   Detailed explanation goes here

%global simDisplayName
global illUnits
virType='crit200';

nboot=1000;
fffun = @(x)(sum(x)./length(x)); % function for bootstrap

radBinEdges=0:0.2:3.6;
mainSequenceType='all';
massThresh=10^9; % threshold for *stellar* mass
%hostBins=[10.7 12:16];

%% parse arguments
i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case {'crit','crit200','200crit'}
            virType='crit200';
        case {'mean','mean200','200mean'}
            virType='mean200';
        case {'500','crit500','500crit'}
            virType='crit500';
        case {'radbinedges','radbins'}
            i=i+1;
            radBinEdges=varargin{i};
        case {'hostbinedges','hostbins'}
            i=i+1;
            radBinEdges=varargin{i};
            %         case 'mask'
            %             i=i+1;
            %             sampleMask=varargin{i};
        case 'nboot'
            i=i+1;
            nboot=varargin{i};
        case 'mainsequence'
            i=i+1;
            mainSequenceType=varargin{i};
        case {'massthresh', 'massthreshold'}
            i=i+1;
            massThresh=varargin{i};
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i});
    end
    i=i+1;
end

%% Perliminaries

subsInfo = illustris.infrastructure.build_sub_fof_connection(subs,fofs);
% define stellar mass
galMass= illustris.utils.get_stellar_mass(subs); % stellar mass within 2*rhalf
% define ssfr
ssfrBase=double(illustris.utils.calc_ssfr(subs,'base',1e-15));

%% define sample
sampleMask=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh);
centrals=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'centrals');
sats=illustris.infrastructure.generateMask('subs',subs','fofs',fofs,'mass',massThresh,'sats');

cenMask=centrals(sampleMask);
satMask=sats(sampleMask);

gmass=galMass(sampleMask);
ssfr=ssfrBase(sampleMask);

%% generate 'main-sequence'
switch lower(mainSequenceType)
    case 'central'
        samMask=cenMask;
        
    case 'all'
        samMask=sampleMask;
end

msMask=samMask & galMass<=10.^10.2;

mainSeq=mk_main_sequence(log10(galMass(msMask)),...
    log10(ssfrBase(msMask)),'bins',8.3:0.15:12.6,...
    'fac',2,'median');

%% extrapolate main sequence relation and define quenched line
msX=log10(massThresh):0.05:12.5;
msY=interp1(mainSeq.msX,mainSeq.msY,msX,'linear','extrap');
qVal=interp1(mainSeq.msX,mainSeq.msY,log10(gmass),'linear','extrap')-1;% quenched is defined as 1 dex below main sequence
qMask=log10(ssfr)<=qVal; % this tells you who is quenched

res.quenchedMask=qMask;
res.mainSeq=mainSeq;
res.sampleMask=sampleMask;

res.qFracGlobal=sum(qMask)./length(qMask);
res.qFracGlobalCentrals=sum(qMask(cenMask))./length(qMask);
res.qFracGlobalSats=sum(qMask(satMask))./length(qMask);


%% define virial parameters
switch lower(virType)
    case 'crit200'
        rvir=double(fofs.Group_R_Crit200(subsInfo.hostFof(sampleMask)+1)); % host rvir for each galaxy
        mvir=double(fofs.Group_M_Crit200(subsInfo.hostFof(sampleMask)+1)).*illUnits.massUnit;  % host mvir for each galaxy
    case 'mean200'
        rvir=double(fofs.Group_R_Mean200(subsInfo.hostFof(sampleMask)+1));
        mvir=double(fofs.Group_M_Mean200(subsInfo.hostFof(sampleMask)+1)).*illUnits.massUnit;
    case 'crit500'
        rvir=double(fofs.Group_R_Crit500(subsInfo.hostFof(sampleMask)+1));
        mvir=double(fofs.Group_M_Crit500(subsInfo.hostFof(sampleMask)+1)).*illUnits.massUnit;
end


%% calculate distances
global LBox
radPosition = findDistance(double(subs.SubhaloPos(:,sampleMask)),double(fofs.GroupPos(:,subsInfo.hostFof(sampleMask)+1)),...
    LBox,3)./rvir;

binIndex=discretize(radPosition,radBinEdges);


%% build Quenched fraction profiles

%% split by hostmass
hmask=false(4,length(mvir));
hmask(1,:)=log10(mvir)>=10.7 & log10(mvir)<12;
hmask(2,:)=log10(mvir)>=12 & log10(mvir)<13;
hmask(3,:)=log10(mvir)>=13 & log10(mvir)<14;
hmask(4,:)=log10(mvir)>=14 & log10(mvir)<15;

res.hostHeader='Host Mass bins: 10.7, 12, 13, 14, 15';


%% split by stellarmass - 4 groups
lgMass=log10(gmass);
smask4=false(4,length(gmass));
smask4(1,:)=lgMass>=8 & lgMass<9;
smask4(2,:)=lgMass>=9 & lgMass<10;
smask4(3,:)=lgMass>=10 & lgMass<11;
smask4(4,:)=lgMass>=11 ;
res.byHostStar4.massHeader='stellar mass bins: 8 9 10 11 and above';

%% split by stellarmass - 2 groups
smask2=false(2,length(gmass));
smask2(1,:)=lgMass>=8 & lgMass<10.5;
smask2(2,:)=lgMass>=10.5;
res.byHostStar2.massHeader='stellar mass bins: 8 10.5 and above';

%% by Host mass only
byHost.qfrac=zeros(length(radBinEdges)-1,4);
byHost.bsci=zeros(length(radBinEdges)-1,2,4);
byHost.xMed=zeros(length(radBinEdges)-1,4);
byHost.count=zeros(length(radBinEdges)-1,4);

for j=1:4  %loop over host mass
    mask1=squeeze(hmask(j,:));
    
    % quenched fraction of centrals
    byHost.qFracCen(j)=sum(qMask(cenMask & mask1))./sum(cenMask & mask1);
    byHost.cenCount(j)=sum(cenMask & mask1);
    
    % quenched fraction profiles of satellites
    mask2=mask1 & satMask;
    for k=1:length(radBinEdges)-1
        mask3=mask2 & binIndex==k;
        
        qf=qMask(mask3);
        byHost.qfrac(k,j)=sum(qf)./length(qf);
        byHost.xMed(k,j)=median(radPosition(mask3));
        byHost.count(k,j)=sum(mask3);
        
        % generate 95% confidence interval from bootstrap
        if sum(mask3)>1
            byHost.bsci(k,:,j)=bootci(nboot,{fffun,qf});
        end
    end
    
end

%% by host Mass and stellar mass - 4 bins

for i=1:4  % in stellar mass bins
    byHostStar4(i).qfrac=zeros(length(radBinEdges)-1,4);
    byHostStar4(i).bsci=zeros(length(radBinEdges)-1,2,4);
    byHostStar4(i).xMed=zeros(length(radBinEdges)-1,4);
    byHostStar4(i).count=zeros(length(radBinEdges)-1,4);
    
    for j=1:4  %loop over host mass
        mask1=squeeze(hmask(j,:)) & squeeze(smask4(i,:));
        
        % quenched fraction of centrals
        byHostStar4(i).qFracCen(j)=sum(qMask(cenMask & mask1))./sum(cenMask & mask1);
        byHostStar4(i).cenCount(j)=sum(cenMask & mask1);
        
        % quenched fraction profiles of satellites
        mask2=mask1 & satMask;
        for k=1:length(radBinEdges)-1
            mask3=mask2 & binIndex==k;
            
            qf=qMask(mask3);
            byHostStar4(i).qfrac(k,j)=sum(qf)./length(qf);
            byHostStar4(i).xMed(k,j)=median(radPosition(mask3));
            byHostStar4(i).count(k,j)=sum(mask3);
            
            % generate 95% confidence interval from bootstrap
            if sum(mask3)>1
                byHostStar4(i).bsci(k,:,j)=bootci(nboot,{fffun,qf});
            end
        end
        
    end
    
end

%% by host Mass and stellar mass - 2 bins

for i=1:2  % in stellar mass bins
    byHostStar2(i).qfrac=zeros(length(radBinEdges)-1,4);
    byHostStar2(i).bsci=zeros(length(radBinEdges)-1,2,4);
    byHostStar2(i).xMed=zeros(length(radBinEdges)-1,4);
    byHostStar2(i).count=zeros(length(radBinEdges)-1,4);
    
    for j=1:4  %loop over host mass
        mask1=squeeze(hmask(j,:)) & squeeze(smask2(i,:));
        
        % quenched fraction of centrals
        byHostStar4(i).qFracCen(j)=sum(qMask(cenMask & mask1))./sum(cenMask & mask1);
        byHostStar4(i).cenCount(j)=sum(cenMask & mask1);
        
        % quenched fraction profiles of satellites
        mask2=mask1 & satMask;
        for k=1:length(radBinEdges)-1
            mask3=mask2 & binIndex==k;
            
            qf=qMask(mask3);
            byHostStar4(i).qfrac(k,j)=sum(qf)./length(qf);
            byHostStar4(i).xMed(k,j)=median(radPosition(mask3));
            byHostStar4(i).count(k,j)=sum(mask3);
            
            % generate 95% confidence interval from bootstrap
            if sum(mask3)>1
                byHostStar4(i).bsci(k,:,j)=bootci(nboot,{fffun,qf});
            end
        end
        
    end
    
end

%%
res.byHost=byHost;
res.byHostStar2=byHostStar2;
res.byHostStar4=byHostStar4;

end
