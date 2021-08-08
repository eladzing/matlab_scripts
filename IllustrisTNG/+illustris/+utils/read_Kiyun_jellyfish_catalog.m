function res=read_Kiyun_jellyfish_catalog(snap,varargin )
%WRITE_CATALOG write data catalog to hdf files in the same structure as the
%other TNG values
%   Detailed explanation goes here

global POSTPROCESSING
%global simName

%% defaults
folderName='flags';
fullPath='';

% if ~exist(snao
% snap='099';
%sim='';

fullName=[POSTPROCESSING '/' folderName '/subhalo_jellyfish_flags_0' num2str(snap) '.hdf5']; 


%% parse arguments
% i=1;
% while(i<=length(varargin))
%     switch lower(varargin{i})
%         case{'path','catpath','catalogpath','pppath'}
%             i=i+1;
%             catPath=varargin{i};
%             if strcmpi(catPath,'default')
%                 catPath=myPOSTPROCESSING;
%             end
%         case{'folder','foldername'}
%             i=i+1;
%             folderName=varargin{i};
%         case{'snap'}
%             i=i+1;
%             snap=varargin{i};
%         case{'simname','sim'}
%             i=i+1;
%             sim=varargin{i};
%         case{'fullname','fullpath'}
%             i=i+1;
%             fullPath=varargin{i};
%         otherwise
%             error('write_catalog - illegal argument: %s',varargin{i})
%     end
%     i=i+1;
% end
% 
% if ~exist('catalogName','var')
%    
%    [fullName,fullPath]=uigetfile({'*.hdf5'},'Select Catalog',[myPOSTPROCESSING '/'],'multiselect','off');
% 
% else
%     
%     
%     if isempty(folderName)
%         folderName=catalogName;
%     end
%     
%     if isempty(fullPath)
%         fullPath=[catPath '/' folderName '/'];
%     end
%     
%     
%     % build correct file name
%     if isempty(sim)
%        
%         sim=simName;
%     end
%     
%     if ~ischar(snap)
%         snap=num2str(snap,'%03d');
%     end
%     
%     
%     fullName=[catalogName '_' snap '.hdf5'];
%     
%     
%     
% end

%% read info on catalog

outputStruct=struct;

infoStruct=h5info(fullName);

len=length(infoStruct.Groups.Datasets);

for i=1:len
    name=infoStruct.Groups.Datasets(i).Name;
    
    arr=h5read(fullName,['/Subhalo/' name]);
    outputStruct.(name)=arr';
end


if len==1
    res=outputStruct.(name);
elseif len>1
    res=outputStruct;
else
    error('read_catalog - strange value ofr mo. of datasets in file: %d',len)
       
end


