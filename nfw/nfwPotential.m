function res = nfwPotential(r,Mv,cv)
%nfwMass  NFW mass profile
%   nfw mass profile needs at least 2 parameters Mvir and cv. r must be in
%   units of Rvir



res=Mv.*nfwA(r,cv)./nfwA(1,cv);

end





end
% i=1;
% while i<=length(varargin)
%     switch lower{varargin{i}}
%         case {'rv','rvir'}
%             i=i+1;
%             Rv=varargin{i};
%         case {'cv','cvir'}
%             i=i+1;
%             cv=varargin{i};
%         case {'rs','rscale'}
%             i=i+1;
%             rs=varargin{i};
%         case {'deltavir'}
%             i=i+1;
%             deltaVir=varargin{i};
%         case {'rho_mean'}
%             i=i+1;
%             rhoMean=varargin{i};
%         case{'zred'}
%             i=i+1;
%             zred=varargin{i};
%         otherwise
%             error('nfwMass - Illegal argument: %s',varargin{i}
%     end
%     i=i+1;
% end
% 
% 

