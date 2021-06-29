function draw_circle_boxx(fn, R, clr)

global hub
%global zred

R=R*hub;

[Sx, Sy , ~] = cylinder(R, 100);

figure(fn);
hold on; 
plot(Sx(:), Sy(:),'Color',clr,'Linewidth',2);
hold off;