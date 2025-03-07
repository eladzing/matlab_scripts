function new_env(hal,varargin) %%typ,aexparg)
% Utility function to set up the FILE_FORMAT and FILE_FORMAT_SPHERE
% environment variables, used globally to locate and load cluster data.
%
% @param halo  Cluster Name
% @param varargin - optiopnal arguments:
%         - sim type: csf / adiabatic
%         - redshift: a1 / a06
%         - OS - Linux / win  - affects path to data files 

%error(nargchk(3,3,nargin,'struct'));

hal=upper(hal);
import cluster
if nargin<1 
    error('new_env: must give cluster id');       
end

clusternames = {'CL101','CL102','CL103','CL104','CL105','CL106','CL107',...
    'CL3','CL5','CL6','CL7','CL9','CL10','CL11','CL14','CL24'};

HALO_BASE_PATH_WIN = 'M:/kravtsov';
%HALO_BASE_PATH_WIN = 'R:/kravtsov/datacubes';
HALO_BASE_PATH_LINUX = '/run/media/zinger/MyPassport/kravtsov';
HALO_BASE_PATH = HALO_BASE_PATH_LINUX;

DEFAULT_PRINTOUT_DIR_WIN='C:/Users/eladzing/Documents/cluster/printout';
DEFAULT_PRINTOUT_DIR_LINUX='/home/zinger/workProjects/cluster/printout'; 
%/home/eladzing/OneDrive/cluster/printout';

if isa(hal,'numeric')
    switch hal
        case{101,102,103,104,105,106,107,3,5,6,7,9,10,11,14,24}
            halo=sprintf('CL%i',hal);
        otherwise
            error('new_env: illegal halo id: %d',hal);
    end
elseif isa(hal,'char')
    switch hal
        case clusternames
            halo=hal;
        otherwise
            error('new_env: illegal halo name: %s',hal);
    end  
else
    error('new_env: illegal halo argument');
end

global zred
global ClusterIndex
global aexpn
global aexp
global hub
global SIMTYPE
global VCM
global DEFAULT_PRINTOUT_DIR
global Omm;
global Oml;

DEFAULT_PRINTOUT_DIR=DEFAULT_PRINTOUT_DIR_LINUX;


global NCELL
NCELL=256;

% DEFAULT values
SIMTYPE='CSF';
aexpn='a1';
aexp=1.0;
c1=true;
c2=true;

i=1;
while i<=length(varargin)
     switch varargin{i}
         case {'csf','cooling'}
             SIMTYPE='CSF';    
             c1=true; 
         case {'adiabatic','ad'}
             SIMTYPE='Adiabatic';
             c1=false;
             %halo=sprintf('%s%s',halo,'a'); WHAT IS THIS SUPPOSED TO DO?
         case {'a1','a0','z0','1','0'}
             aexpn='a1';
             aexp=1.0;
             c2=true;
         case {'a06','z06','a0.6','0.6'}
             aexpn='a06';      
             aexp=0.626;         
             c2=false;
         case {'win','WIN','windos','WINDOWS'}
             HALO_BASE_PATH = HALO_BASE_PATH_WIN;
             DEFAULT_PRINTOUT_DIR=DEFAULT_PRINTOUT_DIR_WIN;
         case {'lin','LIN','linux','LINUX'}
             HALO_BASE_PATH = HALO_BASE_PATH_LINUX;
             DEFAULT_PRINTOUT_DIR=DEFAULT_PRINTOUT_DIR_LINUX;
         case 'path'
             i=i+1;
             HALO_BASE_PATH = varargin{i};
         otherwise
             error('new_env: illegal argument %s',varargin{i})
     end
     i=i+1;
end


clusterdirs = clusternames;

haloformats_adiabatic_a06={'%s_a0.627L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.627L%dMpc.dat',...
    '%s_a0.627L%dMpc.dat','%s_a0.627L%dMpc.dat','%s_a0.627L%dMpc.dat','%s_a0.627L%dMpc.dat',...
    '%s_a0.627L%dMpc.dat','%s_a0.627L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.626L%dMpc.dat',... 
    '%s_a0.627L%dMpc.dat','%s_a0.625L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.625L%dMpc.dat',... 
    '%s_a0.625L%dMpc.dat'};
haloformats_adiabatic_a1={'%s_a1.000L%dMpc.dat','%s_a0.999L%dMpc.dat','%s_a1.000L%dMpc.dat',...
    '%s_a1.000L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a1.000L%dMpc.dat',...
    '%s_a1.000L%dMpc.dat','%s_a0.999L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a1.000L%dMpc.dat',...
    '%s_a0.999L%dMpc.dat','%s_a1.001L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a0.999L%dMpc.dat',...
    '%s_a0.999L%dMpc.dat'};

haloformats_csf_a06={'%s_a0.625L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.625L%dMpc.dat',...
    '%s_a0.626L%dMpc.dat','%s_a0.627L%dMpc.dat','%s_a0.627L%dMpc.dat','%s_a0.626L%dMpc.dat',...
    '%s_a0.626L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.626L%dMpc.dat',... 
    '%s_a0.627L%dMpc.dat','%s_a0.627L%dMpc.dat','%s_a0.626L%dMpc.dat','%s_a0.625L%dMpc.dat',... 
    '%s_a0.626L%dMpc.dat'};
haloformats_csf_a1={'%s_a1.000L%dMpc.dat','%s_a1.001L%dMpc.dat','%s_a1.001L%dMpc.dat',...
    '%s_a0.999L%dMpc.dat','%s_a0.999L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a1.000L%dMpc.dat',...
    '%s_a1.000L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a1.001L%dMpc.dat','%s_a1.000L%dMpc.dat',...
    '%s_a0.999L%dMpc.dat','%s_a1.000L%dMpc.dat','%s_a0.999L%dMpc.dat','%s_a1.001L%dMpc.dat',...
    '%s_a0.999L%dMpc.dat'};

global relaxedList
global unrelaxedList
relaxedList=logical([0 0 0 1 0 0 0 1 1 0 1 0 1 0 1 0]);
unrelaxedList=~relaxedList;


haloidx = find(strcmp(clusternames, halo));
ClusterIndex=haloidx;
if (isempty(haloidx))
    disp(clusternames)
    error('Invalid cluster name!  Must be one of the above');
end

zred=1/aexp-1

if c1
    if c2
        haloformats= haloformats_csf_a1;
    else
        haloformats= haloformats_csf_a06;
    end
else
    if c2
        haloformats= haloformats_adiabatic_a1;
    else
        haloformats= haloformats_adiabatic_a06;
    end
end



global HALO_PATH;
HALO_PATH = sprintf('%s/%s/%s', HALO_BASE_PATH, clusterdirs{haloidx},SIMTYPE);

global FILE_FORMAT;
FILE_FORMAT = sprintf('%s/%s', HALO_PATH, haloformats{haloidx});
global FILE_FORMAT_SPHERE;
FILE_FORMAT_SPHERE = sprintf('%s/%s_%s_%s', HALO_PATH, '%s_sphere',aexpn,'%d.mat');
global FILE_FORMAT_MAT;
FILE_FORMAT_MAT = sprintf('%s/%s_%s_%s', HALO_PATH, '%s',aexpn,'%d.mat');
global PROFILE_FILE;
PROFILE_FILE = sprintf('%s/profile_%s.dat', HALO_PATH,aexpn); 
global FILE_FORMAT_XRAYPROJ;
FILE_FORMAT_XRAYPROJ= sprintf('%s/%s', HALO_PATH, haloformats{haloidx});
global FILE_FORMAT_XRAYCUBE;
FILE_FORMAT_XRAYCUBE= sprintf('%s/%s', HALO_PATH, haloformats{haloidx});

global CLUSTER;
CLUSTER=halo
global RELAXED
RELAXED=relaxedList(haloidx);

% find default vcm
%load(sprintf('matlab/mat_files/vcm_stack_%s.mat',aexpn),'vcm_stack');
load(sprintf('mat_files/vcm_stack_%s.mat',aexpn),'vcm_stack');
hlist=vcm_stack{1};vcmlist=vcm_stack{3};
for i=1:length(hlist)
    if(strcmp(halo,hlist{i}))
        VCM=vcmlist(i,:);
    end
end
clear vcm_stack vcmlist hlist


%POT_TEMP_BASE_PATH='/~/work/sshfs/sungate/poteng/data';
%global POT_FILE_FORMAT
%POT_FILE_FORMAT = sprintf('%s/%s/%s/%s',POT_TEMP_BASE_PATH,CLUSTER,SIMTYPE,haloformats{haloidx})
%global POT_FILE_FORMAT_SPHERE
%POT_FILE_FORMAT_SPHERE = sprintf('%s/%s/%s/%s_%s_%s',POT_TEMP_BASE_PATH,CLUSTER,SIMTYPE,'%s_sphere',aexpn,'%d.mat')


[~, hu, om0, ol0, ~,~,~] = read_header(1);
%[aexp , ~, om0 , ~, ~,~,~,~] = read_header(1);
Omm=om0;
hub=hu;
Oml=ol0;


