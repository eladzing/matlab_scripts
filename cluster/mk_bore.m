
function bore=mk_bore(cube,varargin)%  weight,brindX,brindY,direc)


weight=ones(size(cube));
smoot=2;
i=1;
indFlag=false(1,3);

while i<=length(varargin)
    switch lower(varargin{i})
        case {'weight','wt'}
            i=i+1;
            weight=varargin{i};
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
        case 'zind'
            i=i+1;
            zind=varargin{i};
            indFlag(3)=true;
        otherwise
            error('MK_BORE - Illegal argument; %s',varargin{s});
    end
    i=i+1;
end

cube=cube.*weight;

switch lower(direc)
    case 'x'
        if all(indFlag([2 3]))
            ind1=1:size(cube,1);
            ind2=yind-smoot:1:yind+smoot;
            ind3=zind-smoot:1:zind+smoot;
        else
            error('MK_BORE - index direction mismatch: %s',direc);
        end
        s1=2;
        s2=3;
    case 'y'
        if all(indFlag([1 3]))
            ind2=1:size(cube,2);
            ind1=xind-smoot:1:xind+smoot;
            ind3=zind-smoot:1:zind+smoot;
        else
            error('MK_BORE - index direction mismatch: %s',direc);
        end
        s1=1;
        s2=3;
    case 'z'
        if all(indFlag([1 2]))
            ind3=1:size(cube,3);
            ind2=yind-smoot:1:yind+smoot;
            ind1=xind-smoot:1:xind+smoot;
        else
            error('MK_BORE - index direction mismatch: %s',direc);
        end
        s1=1;
        s2=2;
    otherwise
        error('MK_BORE - Illegal direction: %s',direc);
end

bore=sum(sum(cube(ind1,ind2,ind3),s1),s2)./sum(sum(weight(ind1,ind2,ind3),s1),s2);

% switch direc
%     case {'x','X','1','yz','YZ'}
%         bore=sum(sum(cube(xind,yind,zind),2),3)./sum(sum(weight(xind,yind,zind),2),3);
%     case {'y','Y','2','zx','ZX'}
%         bore=transpose(sum(sum(cube(brindX,:,brindY),1),3)./sum(sum(weight(brindX,:,brindY),1),3));
%     case {'z','Z','3','xy','XY'}
%         bore=squeeze(sum(sum(cube(brindX,brindY,:),1),2)./sum(sum(weight(brindX,brindY,:),1),2));
% end

end
