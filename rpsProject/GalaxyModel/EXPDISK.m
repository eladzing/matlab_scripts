%% defining a Hernquist halo object.
% all radial positions are assumed to be in units of the virial radius

classdef EXPDISK
    properties
        Md
        Rd
        Sigma
        type
        Rhalf
        GG=4.2997e-06; % gravitational constant in Msun kpc km/sec
        header='Exponential disk model in units of solarMass, Kpc, km/sec';
    end
    methods
        % constructor
        function obj=EXPDISK(varargin)
        
        M=[];
        R=[];
        Sig=[];
        type=[];
        Rhalf=[];
        
        i=1;
        while  i<=length(varargin)
            switch(lower(varargin{i}))
                case{'mass','m','md','ms','mg'}
                    i=i+1;
                    M=varargin{i};
                case{'radius','r','rd','rs','rg','rad'}
                    i=i+1;
                    R=varargin{i};
                case{'halfmassradius','rhalf','rh','halflightradius'}
                    i=i+1;
                    Rhalf=varargin{i};
                    
                case{'sig','sigma','sigmad'}
                    i=i+1;
                    Sig=varargin{i};
                case{'type'}
                    i=i+1;
                    type=varargin{i};
                otherwise
                    error('EXPDISK - Illegal argument: %s',varargin{i});
            end
            i=i+1;
        end
        
        
        if isempty(R) && ~isempty(Rhalf)
            R=Rhalf./EXPDISK.effective_rad(0.5);
        end
               
        
        if ~isempty(M)
            if ~isempty(R)
                Sig=M./(2*pi*R.^2);
            elseif ~isempty(Sig)
                R=sqrt(M./(2*pi*Sig)); 
            else 
                error('EXPDISK - something wrong  - not enough arguments given, R or Sigma missing (but M OK)');
            end
        elseif ~isempty(R) && ~isempty(Sig)
            M=Sig.*(2*pi*R.^2);
        else
                error('EXPDISK - something wrong  - not enough arguments given, R or Sigma missing (and no M)');
        end
        
        if isempty(Rhalf)
            Rhalf=EXPDISK.effective_rad(0.5).*R;
        end
        
                
        obj.Md=M;
        obj.Rd=R;
        obj.Sigma=Sig;
        obj.type=type;
        obj.Rhalf=Rhalf;
        
        end
        
        % enclosed mass profile
        function res=mass(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'pc'
                    rr=rr./(obj.Rd*1000);
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
        
                
        res=obj.Md.*EXPDISK.massProf(rr);
                
        end
        
        % surface density profile in solarmass/kpc^2
        function res=density(obj,rr,runit)
        
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'pc'
                    rr=rr./(obj.Rd*1000);
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
        
        res=obj.Sigma.*exp(-rr);
                  
        end
        
        function res=rho(obj,rr,runit)
        res=density(obj,rr,runit);
        end
        
                
        % potential in disk plane 
        function res=potential(obj,rr,runit)
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'pc'
                    rr=rr./(obj.Rd*1000);
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
        
        res=-pi.*obj.GG*obj.Sigma.*obj.Rd.*EXPDISK.dfunc(rr);
        
        end
        
        % circular velocity 
        function res=vcirc(obj,rr,runit)
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'pc'
                    rr=rr./(obj.Rd*1000);
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
       
        res=sqrt(pi*obj.GG.*obj.Sigma.*obj.Rd.*rr.^2.*EXPDISK.bfunc(rr));
        
        end

        
        
        % gravitaional accelaration
        function res=ggrav(obj,rr,runit)
        if exist('runit','var')
            switch(lower(runit))
                case 'kpc'
                    rr=rr./obj.Rd;
                case 'pc'
                    rr=rr./(obj.Rd*1000);
                case 'mpc'
                    rr=rr./(obj.Rd/1000);
            end
        end
       
        res=-1*pi*obj.GG.*obj.Sigma.*rr.*EXPDISK.bfunc(rr);
        
        end

    end
    
    methods (Static)
        
        % function for finding the effective radius for a given mass
        % fraction of an exponential disk in units of rd
        function res=effective_rad(eff)
        
        x=0:1e-4:10;
        y=1-exp(-x).*(1+x);
        
        res=interp1(y,x,eff);
        
        end
        
        % mass profile function in units of Md. rr is in units of Rd. 
        function res=massProf(rr)
        res=1-exp(-rr).*(1+rr);
        end
        
        
        % function to calculate the 'bessel' part of the
        % gravitational acceleration of an exponential disk.
        % r is in units of the scale radius rd 
        function res=bfunc(r)

        res=besseli(0,0.5.*r).*besselk(0,0.5.*r)-besseli(1,0.5.*r).*besselk(1,0.5.*r);

        end
        
        %% auxilary function to calculate the 'bessel' part of the
        % gravitational potential of an exponential disk.
        % r is in units of the scale radius rd
        function res=dfunc(r)
                        
        res=r.*(besseli(0,0.5.*r).*besselk(1,0.5.*r)-besseli(1,0.5.*r).*besselk(0,0.5.*r));

        end
        
        
    end
end
