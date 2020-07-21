function mask = mk_splashback_mask(type,param,flavor,radParam)
%MK_SPLASHBACK_MASK Create mask for identifyng centrals who were Satellites
%not long ago
%   Makes use of existing data structures to create a mask
%   Must decide ihow far back to look to see if it was a satellite
%   Arguments:
%       type - lookback by either by redshift or by time
%       param - how far back (in redshift or time
%       flavor - use radial distance parameter, subhalo identification or
%       both. default is 'both'

% defualts
global simDisplayName
global DEFAULT_MATFILE_DIR
satFlag=true;
radFlag=true;

% check which measures are to be used
if exist('flavor','var')
    switch(lower(flavor))
        case {'subhalo','issat','sat'}
            radFlag=false;
        case {'distance','radius','rad'}
            satFlag=false;
        case {'both'}
        otherwise
            error('MK_SPLASHBACK_MASK - unknown flavor')
    end
end

% controls
if ~ischar(type)
    error('MK_SPLASHBACK_MASK - "type" must by a string')
end


switch(lower(type))
    case{'redhsift','z','zred'}
        ftype='Redshift';
        fname=sprintf('%s/splashbackCatalog_by%s_%s',DEFAULT_MATFILE_DIR,ftype,simDisplayName);
        
        load(fname,'splashbackRedshift');
        
        catalog=splashbackRedshift;
        clear splashbackRedshift
        %zredFlag=true;
        
        pInd=find(catalog.redshift(1:end-1)<param);
        
        fprintf('masking for satellits since z=%s \n',...
            num2str(catalog.redshift(pInd(end)+1)));
        
    case{'time','lookback','t'}
        ftype='Time';
        fname=sprintf('%s/splashbackCatalog_by%s_%s',DEFAULT_MATFILE_DIR,ftype,simDisplayName);
        
        load(fname,'splashbackTime');
        
        catalog=splashbackTime;
        clear splashbackTime
        
        %timeFlag=true;
        
        pInd=find(catalog.Lookback(1:end-1)<param);
        
        fprintf('masking for satellites in the past %s Gyr \n',...
            num2str(catalog.Lookback(pInd(end)+1)));
        
    otherwise
        error('MK_SPLASHBACK_MASK - Illegal option for "type": %s',type)
end

len=length(catalog.isSat);
if satFlag
    satMask=catalog.isSat(pInd,:)==1;
else
    satMask=true(length(pInd),len);
end

if radFlag
    
    % determine radial definition of too far  (in units of R_200,c
    if exist('radParam','var')
        str={'5' '10' '25' '50'};
        ii=find(radParam>=[0.05 0.1 0.25 0.5],1,'last');
        radField=['isFar_' str{ii} 'prc'];
    else
        radField='isFar_10prc';
        
    end
    k=strfind(radField,'prc');
    fprintf('Satellites defined beyond %s %% of R_200,c \n',radField(7:k-1));
        
    
    radMask=catalog.(radField)(pInd,:)==1;
else
    radMask=true(length(pInd),len);
end

mask=any(satMask,1) & any(radMask,1);




end

