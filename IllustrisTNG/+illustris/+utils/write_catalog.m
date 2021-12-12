function  write_catalog(catalog,snap,varargin )
%WRITE_CATALOG write data catalog to hdf files in the same structure as the
%other TNG values
%   Detailed explanation goes here

global  myPOSTPROCESSING

%% defaults
catalogName ='';
folderName='';
fullName='';
catPath=myPOSTPROCESSING;
structFlag=isstruct(catalog);
verboseFlag=false;

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
        case{'name','catname','catalogname'}
            i=i+1;
            catalogName=varargin{i};
        case{'folder','foldername'}
            i=i+1;
            folderName=varargin{i};
        case 'fullname'
            i=i+1;
            fullName=varargin{i};
        case{'verbose','verb','v'}
            verboseFlag=true;
        otherwise
            error('%s - illegal argument: %s',current_function.upper,varargin{i})
    end
    i=i+1;
end

if isempty(catalogName)
    catalogName=inputname(1);
end

if isempty(folderName)
    folderName=catalogName;
end

%% create file if needed

% create directory
fullPath=[catPath '/' folderName];
[s,mess,messid]=mkdir(fullPath);

% build correct file name
%global simName

if isempty(fullName)
    fullName=[catalogName '_' num2str(snap,'%03d') '.hdf5'];
end

%% get info about catalog
if structFlag % Catalog is a structure
    dataSets=fieldnames(catalog);
    %catLength=length(catalog.(dataSets{1}));
else
    catLength=length(catalog);
    dataSets{1}=catalogName;
end



% Printing Data
% Header
% for id = 1:numel(dataSets)
%     dataSet    = char(dataSets(id));
%     what    = {[ char([ catalog.(sprintf([ dataSet ])) ])]};
%     if strcmp(dataSet, 'measure_apertures')
%         what = { catalog.(sprintf([ dataSet ])) };
%     end
%
%     if id == 1
%         h5write(catPath, ['Header/' '' dataSet ''], what);
%     else
%         h5write(catPath, ['Header/' '' dataSet ''], what, 'WriteMode','append');
%     end
% end


%% remove previous versions 

status=move2old([fullPath '/' fullName]);

if ~status
    error('WRITE_CATALOG - something went wrong')
end


%% write catalog 
if structFlag
    for id = 1:numel(dataSets)
        dataSet   = char(dataSets(id));
        catLength=size(catalog.(dataSets{id}));
        
        h5create([fullPath '/' fullName],['/' dataSet],catLength);%,'Datatype',class(catalog.(dataSet)));
        
        h5write([fullPath '/' fullName], ['/' dataSet], catalog.(dataSet));%,'Datatype','single');
        %h5write(catPath, ['' catPs.selection_hType '' '/' '' dataSet ''], sHaloes.(sprintf([ dataSet ])), 'WriteMode', 'append');
        
    end
else
    
    dataSet   = char(dataSets(1));
    h5write([fullPath '/' fullName], ['/' dataSet], catalog, 'WriteMode', 'append');
end


if verboseFlag
    fprintf('*** Finished writing catalog: %s *** \n',[fullPath '/' fullName]);
end

end

