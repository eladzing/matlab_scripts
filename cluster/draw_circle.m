function draw_circle(fn, R)

[Sx Sy Sz] = cylinder(R, 100);

figure(fn);
hold on; 
plot(Sx(:)+128.5, Sy(:)+128.5);
hold off;