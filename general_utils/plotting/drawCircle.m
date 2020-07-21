function drawCircle(figHandle, radius, varargin)
%drawCircle  draw a circle of a given radius in a given figure. 
% color and width can be set;

clr='k';
lw=2;
center=[0 0];
lStyle='-';

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
        case{'clr','color'}
            i=i+1;
            clr=varargin{i};
        case{'width','linewidth','lw'}
            i=i+1;
            lw=varargin{i};
        case 'center'
            i=i+1;
            center=varargin{i};
        case{'linestyle','style','line'}
            i=i+1;
            lStyle=varargin{i};
        otherwise
            error('DRAWCIRCLE - Illegal argument: %s',varargin{i});
            
    end
    i=i+1;
end

[Sx, Sy , ~] = cylinder(radius, 100);

figure(figHandle);
hold on;
plot(Sx(:)+center(1), Sy(:)+center(2),'Color',clr,'Linewidth',lw,'lineStyle',lStyle);
hold off;

end