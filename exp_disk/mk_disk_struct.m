function disc=mk_disk_struct(mstar,rd,varargin)
% function for creating a structure for an galactic disk modeled as an 
% exponentiail disk with an (optional) stellar bulge component. 
% mstar is the total stellar mass
% rd - the stellar component disk scale radius

% the assunmption is that both gas and stars are fit by an exponential disk
% model. there is an optional stellar bulge component
% optional values are:
% fg - gas to stellar mass fraction (default: 0.25); 
% beta - ratio of scale radius of stellar component to gas   (default: 1 )
%
% bulge is assumed to be modeled with a Hernquist porfile
% fb - mass ration of the bulge to disk (stellar only) default=0;
% xi - ratio of the disk scale radius to bulge scale radius default: 0.1

if nargin<2
    error('mk_disk_struct: must enter stellar mass and scale radius')
end

disc.mTotal=mstar;
disc.rd=rd;

% defualt values 
fg= 0.25;
beta=1;
fb=0;
xi=0.1;

i=1;
while i<=length(varargin)
    switch varargin{i}
        case 'fg'
            i=i+1;
            fg=varargin{i};
        case 'beta'
            i=i+1;
            beta=varargin{i};
        case 'fb'
            i=i+1;
            fb=varargin{i};
        case 'xi'
            i=i+1;
            xi=varargin{i};
        otherwise
            error('mk_disk_struct: Illegal argument',varargin{i})
    end
    i=i+1;
end

disc.mstar=disc.mTotal/(1+fb); % find mass in disk stellar component


%% creat observables 
disc.sigma=disc.mstar./(2*pi*disc.rd^2);
disc.sigmaEff=disc.sigma.*sigma_factor('eff');
disc.sigma5090=disc.sigma.*sigma_factor('5090');
disc.BT=fb/(fb+1);

%% create profiles  
r=0.01:0.01:50; % in units of Rd 
disc.rp=r;
rxi=r.*xi; % for use in bulge stuff
% Mass profiles - all are in units of mstar (stellar mass of the disk)
disc.mass.starDisk=exp_disk_mass(r,1);
disc.mass.gasDisk=fg.*exp_disk_mass(r,beta);
disc.mass.totalDisk=disc.mass.starDisk+disc.mass.gasDisk;
disc.mass.bulge=fb.*(rxi/(1+rxi)).^2;  %bulge_mass(r,'fb',fb,'rs',1/xi,'hern');
disc.mass.total=disc.mass.totalDisk+disc.mass.bulge;
disc.mass.totalStars=disc.mass.starDisk+disc.mass.bulge;

% density profiles - all are in units of sigma (mstar/(2 pi Rd^2))
disc.density.starDisk=expdisk_density(r,'stellar');
disc.density.gasDisk=expdisk_density(r,'gas','fg',fg,'beta',beta);
disc.density.totalDisk=disc.density.starDisk+disc.density.gasDisk;
%disc.density.bulge=bulge_density(r,'fb',fb,'rs',1/xi,'hern'); useless in this context 
%disc.desntiy.total=disc.density.totalDisk+disc.density.bulge;

% acceleration profiles - all are in units of G pi Sigma 
disc.accel.starDisk=expdisk_accel(r,'fg',0,'beta',1);
disc.accel.totalDisk=expdisk_accel(r,'fg',fg,'beta',beta);
disc.accel.gasDisk=disc.accel.totalDisk-disc.accel.starDisk;
disc.accel.bulge=2.*fb.*xi^2./(1+rxi).^2;   %bulge_accel(r,'fb',fb,'rs',1/xi,'hern');
disc.accel.total=disc.accel.totalDisk+disc.accel.bulge;

% circular velocity (squared) profiles - all are in units of (G pi Sigma Rd) 
disc.vcSq.starDisk=r.*expdisk_accel(r,'fg',0,'beta',1);
disc.vcSq.totalDisk=r.*expdisk_accel(r,'fg',fg,'beta',beta);
disc.vcSq.gasDisk=disc.vcSq.totalDisk-disc.vcSq.starDisk;
disc.vcSq.bulge=2.*fb.*xi^2.*r./(1+rxi).^2;   %bulge_accel(r,'fb',fb,'rs',1/xi,'hern');
disc.vcSq.total=disc.vcSq.totalDisk+disc.vcSq.bulge;



% potential (squared) profiles - all are in units of mstar 
%disc.mass.starDisk=exp_disk_mass(r,1);
%disc.mass.gasDisk=fg.*exp_disk_mass(r,beta);
%disc.mass.totalDisk=disc.mass.starDisk+disc.mass.gasDisk;
%disc.mass.bulge=bulge_mass(r,'fb',fb,'rs',1/xi,'hern');
%disc.mass.total=disc.mass.totalDisk+disc.mass.bulge;








