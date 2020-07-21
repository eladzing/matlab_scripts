%% load data
global simDisplayName


snap=99; 

global DEFAULT_MATFILE_DIR
load([DEFAULT_MATFILE_DIR '/gasCellSize_snp' num2str(snap) '_' simDisplayName '.mat'])


%% determine mask 

msk=cellSizeStruct.galMask; 
cenMask=cellSizeStruct.galMaskCentrals; 
satMask=cellSizeStruct.galMaskSats;
bkspMask=cellSizeStruct.galMaskBksp;

cenMask2=cenMask & ~bkspMask; 



%% find Mean and standard deviation 


for kk=1:5
    
    switch kk
        case 1
            fld='inGal';
        case 2
            fld='inCGM';
        case 3
            fld='inOut';
        case 4
            fld='inSub';
        case 5
            fld='inAll';
    end
    
    for jj=1:3
        switch jj
            case 1
                tag='All';
            case 2
                tag='SF';
            case 3
                tag='NSF';
        end
        
        % all 
        pnts=cellSizeStruct.(fld).(tag).clMean(msk);
        pnts2=cellSizeStruct.(fld).(tag).clMean2(msk);
        wts=cellSizeStruct.(fld).(tag).clCount(msk);
        wts(pnts==0)=0;
        
        
        meanCellSize.(fld).(tag).meanAll=sum(pnts.*wts)/sum(wts);
        meanCellSize.(fld).(tag).stdAll=sqrt(sum(pnts2.*wts)/sum(wts)-...
            (sum(pnts.*wts)/sum(wts))^2);
        
        % centrals 
        pnts=cellSizeStruct.(fld).(tag).clMean(cenMask);
        pnts2=cellSizeStruct.(fld).(tag).clMean2(cenMask);
        wts=cellSizeStruct.(fld).(tag).clCount(cenMask);
        wts(pnts==0)=0;
        
        meanCellSize.(fld).(tag).meanCen=sum(pnts.*wts)/sum(wts);
        meanCellSize.(fld).(tag).stdCen=sqrt(sum(pnts2.*wts)/sum(wts)-...
            (sum(pnts.*wts)/sum(wts))^2);
        
        % centrals (no back-splash)
        pnts=cellSizeStruct.(fld).(tag).clMean(cenMask2);
        pnts2=cellSizeStruct.(fld).(tag).clMean2(cenMask2);
        wts=cellSizeStruct.(fld).(tag).clCount(cenMask2);
        wts(pnts==0)=0;
        
        meanCellSize.(fld).(tag).meanCen2=sum(pnts.*wts)/sum(wts);
        meanCellSize.(fld).(tag).stdCen2=sqrt(sum(pnts2.*wts)/sum(wts)-...
            (sum(pnts.*wts)/sum(wts))^2);
        
        % Satellites 
        pnts=cellSizeStruct.(fld).(tag).clMean(satMask);
        pnts2=cellSizeStruct.(fld).(tag).clMean2(satMask);
        wts=cellSizeStruct.(fld).(tag).clCount(satMask);
        wts(pnts==0)=0;
        
        meanCellSize.(fld).(tag).meanSat=sum(pnts.*wts)/sum(wts);
        meanCellSize.(fld).(tag).stdSat=sqrt(sum(pnts2.*wts)/sum(wts)-...
            (sum(pnts.*wts)/sum(wts))^2);
        
        % BackSplash
        pnts=cellSizeStruct.(fld).(tag).clMean(bkspMask);
        pnts2=cellSizeStruct.(fld).(tag).clMean2(bkspMask);
        wts=cellSizeStruct.(fld).(tag).clCount(bkspMask);
        wts(pnts==0)=0;
        
        meanCellSize.(fld).(tag).meanBksp=sum(pnts.*wts)/sum(wts);
        meanCellSize.(fld).(tag).stdBksp=sqrt(sum(pnts2.*wts)/sum(wts)-...
            (sum(pnts.*wts)/sum(wts))^2);
        
        
    end
end


