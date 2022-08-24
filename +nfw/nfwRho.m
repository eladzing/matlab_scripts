function res = nfwRho(r,cv)
%nfwMass  NFW mass profile 
%   normalized density answer should be multiplied by Mv/(4*pi/3*Rv^3)
%   r is in units of RV 

res=1./(3.*nfw.nfwA(ones(size(r)),cv).*r.*(1./cv+r).^2);




% i=1;
% while i<=length(varargin)
%     switch lower{varargin{i}}
%         case {'rv','rvir'}
%             i=i+1;
%             Rv=varargin{i};
% case {'cv','cvir'}
%             i=i+1;
%             cv=varargin{i};
% case {'rs','rscale'}
%             i=i+1;
%             rs=varargin{i};
% case {'deltavir'}
%             i=i+1;
%             deltaVir=varargin{i};
% case {'rho_mean'}
%             i=i+1;
%             rhoMean=varargin{i};
%         case{'zred'}
%             zred=
%         
%         otherwise
%             error('nfwMass - Illegal argument: %s',varargin{i}
%     end
%     i=i+1;
% end



end

