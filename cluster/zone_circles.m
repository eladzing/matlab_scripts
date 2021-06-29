function zone_circles(zone) %%zone is a vector of zone radii

%% draws a circle demarking the zones at radii given by the variable 'zone'
% @param zone is in units of R_vir

global hub;
rv=get_rvir().*hub;
for i=1:length(zone)
    if zone(i)==1
         draw_circle_boxx(gcf,rv,'white');       
    else
        draw_circle_boxx(gcf,zone(i)*rv,'black');
    end
end

end

