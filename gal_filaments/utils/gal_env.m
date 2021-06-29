function gal_env(galName,varargin) %%typ,aexparg)
% Utility function to set up the FILE_FORMAT and PATHs 
% environment variables, used globally to locate and load cluster data.
%
% @param galName name of galaxy 
% @param varargin - optiopnal arguments:

   %error(nargchk(3,3,nargin,'struct'));


global DEFAULT_PRINTOUT_DIR
global GAL_NAME
global GAL_BASE_PATH
global GAL_PATH 
global NCELL
global boxSize
global zred
global RVIR

GAL_BASE_PATH = 'M:/nirm';
GAL_NAME = galName;
DEFAULT_PRINTOUT_DIR='C:/Users/eladzing/Documents/filaments/gal_filaments/printout';

GAL_PATH=[ GAL_BASE_PATH '/' GAL_NAME ];

%NCELL=300;




 
