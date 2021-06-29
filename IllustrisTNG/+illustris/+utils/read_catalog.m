function res=read_catalog(catalogName,varargin )
%WRITE_CATALOG write data catalog to hdf files in the same structure as the
%other TNG values
%   Detailed explanation goes here

global POSTPROCESSING
global simName

%% defaults
folderName='';
fullPath='';
snap='099';
sim='';
catPath=POSTPROCESSING;

%% parse arguments
i=1;
while(i<=length(varargin))
    switch lower(varargin{i})
        case{'path','catpath','catalogpath','pppath'}
            i=i+1;
            catPath=varargin{i};
            if strcmpi(catPath,'default')
                catPath=myPOSTPROCESSING;
            end
        case{'folder','foldername'}
            i=i+1;
            folderName=varargin{i};
        case{'snap'}
            i=i+1;
            snap=varargin{i};
        case{'simname','sim'}
            i=i+1;
            sim=varargin{i};
        case{'fullpath'}
            i=i+1;
            fullPath=varargin{i};
        case{'fullname'}
            i=i+1;
            
        otherwise
            error('write_catalog - illegal argument: %s',varargin{i})
    end
    i=i+1;
end

if ~exist('catalogName','var')
   
   [fullName,fullPath]=uigetfile({'*.hdf5'},'Select Catalog',[POSTPROCESSING '/'],'multiselect','off');

elseif ~isempty(fullPath)
    
    fullName=[catalogName '.hdf5'];
    
else
    
    
    if isempty(folderName)
        folderName=catalogName;
    end
    
    %if isempty(fullPath)
    fullPath=[catPath '/' folderName ];
    %end
    
    
    % build correct file name
    if isempty(sim)
       
        sim=simName;
    end
    
    if ~ischar(snap)
        snap=num2str(snap,'%03d');
    end
    
    
%     fullName=[catalogName '_' sim snap '.hdf5'];
fullName=[catalogName '_' snap '.hdf5'];
    
    
    
end

%% read info on catalog

outputStruct=struct;

infoStruct=h5info([fullPath '/' fullName]);
len=length(infoStruct.Datasets);

for i=1:len
    name=infoStruct.Datasets(i).Name;
    
    outputStruct.(name)=h5read([fullPath '/'  fullName],['/' name]);
end


if len==1
    res=outptutStruct.(name);
elseif len>1
    res=outputStruct;
else
    error('read_catalog - strange value ofr mo. of datasets in file: %d',len)
       
end


