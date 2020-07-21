function h=my_subplot(row,col,n)
%subplot('Position',[left bottom width height])
start_left=0.1;
start_bottom=0.1;
width=(1-2*start_left)/col;
height=(1-2*start_bottom)/row;
h=subplot('Position',[(start_left+mod((n-1),col)*width) (start_bottom+floor((n-1)/col)*height) width height]);

set(h,'xtick',[],'ytick',[]);
%axis off;
hold on
if (mod((n-1),col)==0 ) %
%Reset the bottom subplot to have xticks
set(h,'ytickMode', 'auto');
end

if ( floor((n-1)/col)==0 )
    set(h,'xtickMode', 'auto');
end
box on
hold on
end