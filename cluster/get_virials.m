function result=get_virials(varargin)
%% get rvir/r200/r500 from list.
%optional argument can be empty or 200 or 500
% returns an array of 5 values: rvir mvir vvir tvir delta_vir
    
  global CLUSTER 
  global aexpn
  load(sprintf('mat_files/virials_%s',aexpn));

  
  haloidx = strcmp(cluster_virials(:,1), CLUSTER);
    
    
%  Don't want more than 3 printing arugments altogether
numvarargs=length(varargin);
if numvarargs>1
    error('get_rvir: too many printing arguments')
end
    
%set defualt printing flag, printing tag, and output dir
defvals={1};

%assign the optional values
defvals(1:numvarargs)=varargin;

% transfer to easy to use varaibles
[deltaflag]=defvals{:};

switch deltaflag
    case 1
        index=2;
    case 200
        index=3;
    case 500
        index=4;
    otherwise
        error('get_rvir: illegal option for argument')
end

result=cluster_virials{haloidx,index};

%result=virs(1);

end

    
    
    
    