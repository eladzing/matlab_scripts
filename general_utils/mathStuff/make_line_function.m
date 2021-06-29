function [a,b ] = make_line_function( p1,p2 )
% MAKE_LINE_FUNCTION -given two points p1,p2, derive the parameters for the 
% straight line connecting the two: y=a*x+b

a=(p2(2)-p1(2))./(p2(1)-p1(1));
b=(p1(2)*p2(1)-p2(2)*p1(1))./(p2(1)-p1(1));


end

