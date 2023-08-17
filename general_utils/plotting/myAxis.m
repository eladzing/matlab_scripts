function ax = myAxis(varargin)
%MYAXIS default axis configuration
%   Detailed explanation goes here

ax=gca;
%axPos=get(ax.'position');
fontsize=14;

i=1;
while i<=length(varargin)
    switch lower(varargin{i})
        
        case{'font','size','fontsize'}
            i=i+1;
            fontsize=varargin{i};
        
%         case{'pos','position'}
%             i=i+1;
%             axPos=varargin{i};
            
        otherwise
            error('%s - Illegal argument: %s',current_function().upper,varargin{i})
    end
    i=i+1;
end
   
set(ax,'fontsize',fontsize,'TickLabelInterpreter','latex')

end

