
global simDisplayName

if ~skip2write
    %global illUnits
    global cosmoStruct
    
    if readFlag
        
        global DEFAULT_MATFILE_DIR
        load([DEFAULT_MATFILE_DIR '/central_history_splashback_z0_' simDisplayName '.mat'])
        
        
    end
    
    fprintf(' *** Finished reading  *** \n');
    
    
    radThresh=[0.05 0.10 0.25 0.5];
    
    ind=find(done);
    
    %timeFrame=[1 2 5 8 12]
    %   1Gyr    2Gyr  5.2Gyr  8  10.5  12.25
    zz=[0 0.1 0.2  0.5 1 2 4];
    time=[0 1 2 5 10 12];
    tim=redshift2time(zz,'cosmo',cosmoStruct);
    lookback=tim.lookback;
    
    isSatZ=zeros(length(zz)-1,length(done));
    isSatT=zeros(length(time)-1,length(done));
    
    isFarZ=zeros(length(radThresh),length(zz)-1,length(done));
    isFarT=zeros(length(radThresh),length(time)-1,length(done));
    
    isSatZ(:,~done)=-1;
    isFarZ(:,:,~done,:)=-1;
    isSatT(:,~done)=-1;
    isFarT(:,:,~done,:)=-1;
    
    step=10;
    thresh=step;
    for i=1:length(ind)
        
        if floor(100*i/length(ind))>thresh
            fprintf('%s %% completed \n',num2str(thresh));
            thresh=thresh+step;
        end
        
        
        rad=centralHist(ind(i)).radiusToHost./centralHist(ind(i)).hostR200;
        
        isSat0=~centralHist(ind(i)).isCentral';
        
        % find 3 concurrent instances of being not central
        isSatm1=[false isSat0(1:end-1)];
        isSatp1=[isSat0(2:end) false];
        isSatHist=isSat0 & isSatm1 & isSatp1;
        
        zSat=centralHist(ind(i)).zred(isSatHist);
        tim=redshift2time(zSat,'cosmo',cosmoStruct);
        lookBSat=tim.lookback;
        
        
        for k=1:length(zz)-1
            isSatZ(k,ind(i))=any(zSat>zz(k) & zSat<zz(k+1));
        end
        
        for k=1:length(time)-1
            isSatT(k,ind(i))=any(lookBSat>time(k) & lookBSat<time(k+1));
        end
        
        
        % find 3 concurrent instances of being not central
        for j=1:length(radThresh)
            isFar0=rad>radThresh(j);
            isFarm1=[false isFar0(1:end-1)];
            isFarp1=[isFar0(2:end) false];
            isFarHist=isFar0 & isFarm1 & isFarp1;
            
            zFar=centralHist(ind(i)).zred(isFarHist);
            tim=redshift2time(zSat,'cosmo',cosmoStruct);
            lookBFar=tim.lookback;
            
            
            for k=1:length(zz)-1
                isFarZ(j,k,ind(i))=any(zFar>zz(k) & zFar<zz(k+1));
            end
            
            for k=1:length(time)-1
                isFarT(j,k,ind(i))=any(lookBFar>time(k) & lookBFar<time(k+1));
            end
        end
        
    end
end
fprintf(' *** Finished preparing catalog. Now saving  *** \n');

% zz=[0 0.1 0.2  0.5 1 2 4];
% time=[0 1 2 5 10 12];
splashbackRedshift.isSat=isSatZ;
splashbackRedshift.redshift=zz;

%radThresh=[0.05 0.10 0.25 0.5 1];
splashbackRedshift.isFar_5prc=squeeze(isFarZ(1,:,:));
splashbackRedshift.isFar_10prc=squeeze(isFarZ(2,:,:));
splashbackRedshift.isFar_25prc=squeeze(isFarZ(3,:,:));
splashbackRedshift.isFar_50prc=squeeze(isFarZ(4,:,:));

% time=[0 1 2 5 10 12];
splashbackTime.isSat=isSatT;
splashbackTime.Lookback=time;

splashbackTime.isFar_5prc=squeeze(isFarT(1,:,:));
splashbackTime.isFar_10prc=squeeze(isFarT(2,:,:));
splashbackTime.isFar_25prc=squeeze(isFarT(3,:,:));
splashbackTime.isFar_50prc=squeeze(isFarT(4,:,:));





% splashback.isSatT=isSatT;
% splashback.isFarZ=isFarZ;
% splashback.isFarT=isFarT;
% splashback.lookback=time;
% splashback.zredBins=zz;
% splashback.done=done;
% splashback.radThresh=radThresh;

fname=sprintf('splashbackCatalog_byRedshift_%s',simDisplayName);

save([DEFAULT_MATFILE_DIR '/' fname],'splashbackRedshift','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);


fname=sprintf('splashbackCatalog_byTime_%s',simDisplayName);

save([DEFAULT_MATFILE_DIR '/' fname],'splashbackTime','-v7.3')

fprintf(' *** Result saved to: %s *** \n',[DEFAULT_MATFILE_DIR '/' fname]);





