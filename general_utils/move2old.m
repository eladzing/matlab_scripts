function status =move2old(fullName,remove )
%MOVE2OLD Check to see if a file alreadt exists and if so move to 'old'
%subfolder or remove completely if stated.
%   Detailed explanation goes here

if ~exist('remove','var')
    rmFlag=false;
else
    switch lower(remove)
        case{'rm','del','remove','delete','erase'}
            rmFlag=true;
        otherwise
            rmFlag=false;
    end
end
          


%% decompose path and filename 
mark=strfind(fullName,'/');

if isempty(mark)
    pathe='./';
    fileName=fullName;
else
    pathe=fullName(1:max(mark));
    fileName=fullName(max(mark)+1:end);
end
    

%% check to see if file exists

fileThere=exist(fullName,'file');

if fileThere==0
    fprintf('%s - The folder %s does not contain the a file by the name %s. A O K. \n',current_function().upper,pathe,fileName);
    status=true;
elseif fileThere==2
    if rmFlag
        delete(fullName);
        status=true;
    else
        % check to see if 'old' folder exists/ 
        if exist([pathe 'old'],'dir')==0
            mkdir([pathe 'old'])
        end
        
        status=movefile(fullName,[pathe 'old/' fileName]);
    end
else
    error('%s- strange error',current_function().upper);
end
        
        

end

