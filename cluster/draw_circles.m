function draw_circles(zone,markColor) %%zone is a vector of zone radii

%% draws a circle demarking the zones at radii given by the variable 'zone'
% @param zone is in units of Mpc 

%global hub;
%rv=get_rvir().*hub;
for i=1:length(zone)
    draw_circle_boxx(gcf,zone(i),markColor);
 end

end

