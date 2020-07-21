function tff = recreate_TFF(rpos,r200c,a,b,c)
%RECREATE_TFF recreate the free-fall time profile based on the polynomial
%fit
%   must enter the position in halo, the r200,c and the three fit
%   coefficents



rr=log10(rpos./r200c);

tff=10.^(a.*rr.^2+b.*rr+c);



end

