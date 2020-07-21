MPSec = 8;
[r_mesh phi_mesh theta_mesh] = sphere_grid(MPSec);

fff = full_ff(:,1:10:end,1:10:end);
fff1 = reshape(fff, [256, 26^2]);
plot(fff1)
semilogy(fff1)
plot(fff1)
fff = full_ff(:,1:20:end,1:20:end);
fff1 = reshape(fff, [256, 13^2]);
cmaskkk1 = reshape(full_cold_mask(:,1:20:end,1:20:end), [256, 13^2]);
fff2 = fff1;
fff2(~cmaskkk1) = 0;
plot(fff2)
plot(fff2)

SS = S_sphere(4);
SS1 = reshape(SS(:,1:20:end,1:20:end), [256, 13^2]);

RR = squeeze(r_mesh(:,1,1));
loglog(RR/2.2, SS1)
loglog(RR/2.2, sum(SS1,2))

xlim([0 0.15])
ylim([0 0.1])

RRR = RHOG_sphere(8);
RR1 = reshape(RRR(:,1:20:end,1:20:end), [256, 13^2]);
plot(RR/2.2, RR1)