
function bore=mk_bore_xrayProj(proj,varargin)%  weight,brindX,brindY,direc)


% weight=ones(size(proj));
smoot=2;
i=1;
indFlag=false(1,2);

while i<=length(varargin)
    switch lower(varargin{i})
%         case {'weight','wt'}
%             i=i+1;
%             weight=varargin{i};
         case {'dir','direction'}
            i=i+1;
            direc=varargin{i};
        case {'smooth'}
            i=i+1;
            smoot=varargin{i};
        case 'xind'
            i=i+1;
            xind=varargin{i};
            indFlag(1)=true;
        case 'yind'
            i=i+1;
            yind=varargin{i};
            indFlag(2)=true;
%         case 'zind'
%             i=i+1;
%             zind=varargin{i};
%             indFlag(3)=true;
        otherwise
            error('MK_BORE - Illegal argument; %s',varargin{s});
    end
    i=i+1;
end

% proj=proj.*weight;

switch lower(direc)
    case 'x'
        if indFlag(2)
            ind1=1:size(proj,1);
            ind2=yind-smoot:1:yind+smoot;
%             ind3=zind-smoot:1:zind+smoot;
        else
            error('MK_BORE - index direction mismatch: %s',direc);
        end
        s1=1;
%         s2=3;
    case 'y'
        if indFlag(1)
            ind2=1:size(proj,2);
            ind1=xind-smoot:1:xind+smoot;
            %ind3=zind-smoot:1:zind+smoot;
        else
            error('MK_BORE - index direction mismatch: %s',direc);
        end
        s1=2;
%         s2=3;
%     case 'z'
%         if all(indFlag([1 2]))
%             ind3=1:size(proj,3);
%             ind2=yind-smoot:1:yind+smoot;
%             ind1=xind-smoot:1:xind+smoot;
%         else
%             error('MK_BORE - index direction mismatch: %s',direc);
%         end
%         s1=1;
%         s2=2;
    otherwise
        error('MK_BORE_XRAYPROJ - Illegal direction: %s',direc);
end

bore=squeeze(sum(proj(ind2,ind1,:),s1))./(2*smoot+1);

%bore=sum(sum(proj(ind1,ind2,ind3),s1),s2)./sum(sum(weight(ind1,ind2,ind3),s1),s2);

% switch direc
%     case {'x','X','1','yz','YZ'}
%         bore=sum(sum(cube(xind,yind,zind),2),3)./sum(sum(weight(xind,yind,zind),2),3);
%     case {'y','Y','2','zx','ZX'}
%         bore=transpose(sum(sum(cube(brindX,:,brindY),1),3)./sum(sum(weight(brindX,:,brindY),1),3));
%     case {'z','Z','3','xy','XY'}
%         bore=squeeze(sum(sum(cube(brindX,brindY,:),1),2)./sum(sum(weight(brindX,brindY,:),1),2));
% end

end
