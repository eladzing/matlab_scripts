%% defining an NFW halo object.
% all radial positions are assumed to be in units of the virial radius

classdef NFW
    properties
        Mvir
        cvir
        Rvir
        Tvir
        Vvir
        fg
        GG=4.2997e-06; % gravitational constant in Msun kpc km/sec
        header='NFW model in units of solarMass, Kpc, km/sec, refernce density: %s';
    end
    methods
        % constructor
        function obj=NFW(varargin)
        
            global cosmoStruct
            
        fg=0.14;
        mv=[];
        cc=[];
        rv=[];
        cosmoSt=[];
        
        delv=200;
        rhoType='crit';
        zred=0;
        
        i=1;
        while i<=length(varargin)
            switch(lower(varargin{i}))
                case{'fg','fgas'}
                    i=i+1;
                    fg=varargin{i};
                case{'mv','mvir','m'}
                    i=i+1;
                    mv=varargin{i};
                case{'cc','cvir','c','cv'}
                    i=i+1;
                    cc=varargin{i};
                    
                case{'rv','rvir','r'}
                    i=i+1;
                    rv=varargin{i};
                    
                    
                case {'delv','delta','delta_vir','deltavir'}
                    i=i+1;
                    delv=varargin{i};
                    
                    
                case {'mean','rho_mean','rhomean'}
                    rhoType='mean';
                case {'crit','rho_crit','rhocrit'}
                    rhoType='crit';
                    
                case {'cosmostruct','cosmo'}
                    i=i+1;
                    cosmoSt=varargin{i};
                    if ~isstruct(cosmoSt)
                        error('NFW: cosmology structure must be a structure')
                    end
                    
                otherwise
                    error('NFW: illegal argument: %s',varargin{i})
            end
            i=i+1;
        end
        
        if isempty(mv) 
            error('NFW: must enter virial mass');
        end
        
       
        
        if isempty(cosmoSt)
            
            if ~isempty(cosmoStruct)
                
                cosmoSt=cosmoStruct;
            else
                cosmoSt=set_LCDM_cosmology('noshow');
                fprintf('cosmology set to LCDM (Planck2015/2016) \n');
            end
            
            
        end
        
        if isempty(cc)
            cc=cvir_Mvir(mv,0,'hub',cosmoSt.hub);
            fprintf('concentration set by cvir-mvir relation - not random \n') ;
        end
        
        
        if isempty(rv)
            [rv, ~, tv, vv]=calculate_virials('mvir',mv,...
                'cosmoStruct',cosmoSt,'delv',delv','zred',zred,rhoType);
            
        else
            [~, ~, tv, vv]=calculate_virials('mvir',mv,...
                'cosmoStruct',cosmoSt,'delv',delv','zred',zred,rhoType);
        end
        
        obj.header=sprintf(obj.header,rhoType);
        
        obj.Mvir=mv;  %in Msun
        obj.Rvir=1000.*rv; %in kpc
        obj.Tvir=tv; %in kelvin
        obj.Vvir=vv; % in km/sec
        obj.cvir=cc; % unitless
        obj.fg=fg;
        end
        
        
        % enclosed mass profile
        function res=mass(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
        end
        
        res=obj.Mvir.*NFW.Afunc(rr,obj.cvir)./NFW.Afunc(1,obj.cvir);
        
        end
        
        % density profile
        function res=density(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
        end
        
        fac=obj.Mvir/(4*pi*obj.Rvir^3)/NFW.Afunc(1,obj.cvir);
        res=fac./rr./(obj.cvir^(-1)+rr).^2;
        end
        
        function res=rho(obj,rr,runit)
            
            if exist('runit','var')
                
                res=density(obj,rr,runit);
            else
                res=density(obj,rr);
            end
        end
        
        % gas density profile
        function res=gasDensity(obj,rr,runit)
        
        
        res=obj.fg.*density(obj,rr,runit);
        end
        
        
        
        % potential
        function res=potential(obj,rr,runit)
        
         if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
         end
        
             
        
        fac=-1.*obj.GG.*obj.Mvir./obj.Rvir./NFW.Afunc(1,obj.cvir);
        res=fac.*log(1+obj.cvir.*rr)./rr;
        end
        
         % vcirc 
        function res=vcirc(obj,rr,runit)
        
         if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
         end
         
      
                      
        Mr=obj.Mvir.*NFW.Afunc(rr,obj.cvir)./NFW.Afunc(1,obj.cvir);
         
        
        
        res=sqrt(obj.GG.*Mr./(obj.Rvir.*rr));
                
        end
      
        
           % gravitational acceleration
        function res=ggrav(obj,rr,runit)
        
         if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
         end
         
      
                      
        Mr=obj.Mvir.*NFW.Afunc(rr,obj.cvir)./NFW.Afunc(1,obj.cvir);
         
        
        
        res=-1.*obj.GG.*Mr./(obj.Rvir.*rr).^2;
                
        end
      
        
        % escape velocity 
        function res=vesc(obj,rr,runit)
        
         if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
         end
        
        
        res=sqrt(-2.*potential(obj,rr));
        end
        
         % free fall time  
        function res=freefallTime(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
        end
         
        Gyr=3.1558e+16;
        kmkpc=1e5./3.0856e+21;
        
        Mr=obj.Mvir.*NFW.Afunc(rr,obj.cvir)./NFW.Afunc(1,obj.cvir);
        rho=Mr./(4*pi/3.*(rr.*obj.Rvir).^3);   % in solarmass kpc^-3
        
        res=sqrt(3*pi./(32*obj.GG.*rho.*kmkpc^2))./Gyr;
        
        end
        
        % mean density 
        function res=meanDensity(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rvir;
                case 'mpc'
                    rr=rr./(obj.Rvir/1000);
            end
        end
         
        
        Mr=obj.Mvir.*NFW.Afunc(rr,obj.cvir)./NFW.Afunc(1,obj.cvir);
        res=Mr./(4*pi/3.*(rr.*obj.Rvir).^3);   % in solarmass kpc^-3
                
        
        end
        
        
    end
    
    methods (Static)
          % auxillary function for mass profile
        function res=Afunc(r,c)
        
        xx=c.*r;
        
        res=log(1+xx)-xx./(1+xx);
        end
    end
end
