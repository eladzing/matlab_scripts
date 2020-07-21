function reslim=find_inner_reslim(varargin)
%% Function to find the inner limit for plotting based on the minimal resolution
% @param varargin is an optional input :  { boxsize ; units ; no. of cells}   
% function finds the smallest  cell size in proper Mpc of the given boxsize
% If units is given it will recast the result in the units specified 
% this is useful mainly for plotting in units of rvir
% the no. of cells designates the number of cells over which we deem our
% limit to be 
%Don't want more than 3 arugments altogether
boxx=1;

numvarargs=length(varargin);
if numvarargs>3
    error('find_inner_reslim: too many inputs','Only 3 are allowed');
end

%set defualt  fontsize
% { boxsize ; units ; no. of cells}   
defvals={1.0 5};

%assign the optional values
defvals(1:numvarargs)=varargin;

% transfer to easy to use varaibles
[units num]=defvals{:};
global hub
global NCELL
global zred

reslim=boxx/(1+zred)/hub/NCELL/units*num;






