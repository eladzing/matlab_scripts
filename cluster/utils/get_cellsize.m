function cellSize = get_cellsize(boxx,varargin)
%% get_cellsize - calculate cellsize of a given box 
% boxx is in Mpc/h;

global hub ;
global NCELL; 

units; 
cellSizeBase=boxx./NCELL; % in Mpc/h 

if isempty(varargin)
    error('get_cellsize: must enter units')
end
i=1;
while i<=length(varargin)
    switch lower(varargin{i})
          case 'mpc'
        cellSize=cellSizeBase./hub;
        case 'mpc/h'
            cellSize=cellSizeBase;
    case 'cm'
       cellSize=cellSizeBase./hub.*Units.Mpc;
       case 'kpc'
       cellSize=cellSizeBase./hub.*1000;

       case 'kpc/h'
       cellSize=cellSizeBase.*1000;
        
        otherwise 
            error('get_cellsize - Illegal argument: %s',varargin{i});
    end
    i=i+1;
end


end





