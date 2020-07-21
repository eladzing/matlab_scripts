units;
host=NFW('mv',1e14,'cc',5);

IC.x=host.Rvir;
IC.y=0;
IC.vx=0;
IC.vy=0.5.*host.Vvir;


tmax=2*(2*pi*host.Rvir/host.Vvir);

dt=1*tmax/4/1e4;
dtmin=dt/1e3;

orb=orbits.rk4_orbitIntegration(IC,dt,tmax,dtmin,@orbits.rhs_nfw,host);

timeunit=(Units.kpc/Units.km)/(Units.Gyr);

%% plot figures
tim=orb.t.*timeunit;
%% plot orbit 
 figure
 plot(orb.x./host.Rvir,orb.y./host.Rvir,'-')
 grid
 axis equal
 xlabelmine('$X/R_{\mathrm{vir}}$')
 ylabelmine('$Y/R_{\mathrm{vir}}$')
 
 %% plot radius
 figure
 plot(tim,orb.rad./host.Rvir,'.-')
 grid
 xlabelmine('$t\,[\mathrm{Gyr}]$')
 ylabelmine('$r/R_{\mathrm{vir}}$')
 
 %% plot x,y 
 figure
 h(1)=plot(tim,orb.x./host.Rvir,'b.-',...
     'DisplayName','$X$');
 hold on 
 h(2)=plot(tim,orb.y./host.Rvir,'r.-',...
     'DisplayName','$Y$');
 grid
 xlabelmine('$t\,[\mathrm{Gyr}]$')
 ylabelmine('$X/R_{\mathrm{vir}},\,Y/R_{\mathrm{vir}}$')
 hl=legend(h);
 set(hl,'Interpreter','latex','fontsize',14,'location','SouthWest');
 
 %% plot vx,vy 
 figure
 h(1)=plot(tim,orb.vx./host.Vvir,'b.-',...
     'DisplayName','$v_\mathrm{x}$');
 hold on 
 h(2)=plot(tim,orb.vy./host.Vvir,'r.-',...
     'DisplayName','$v_\mathrm{y}$');
 grid
 xlabelmine('$t\,[\mathrm{Gyr}]$')
 ylabelmine('$v_\mathrm{x}/V_{\mathrm{vir}},\,v_\mathrm{y}/V_{\mathrm{vir}}$')
 hl=legend(h);
 set(hl,'Interpreter','latex','fontsize',14,'location','SouthWest');
 
 %% plot phase space 
 figure
 h(1)=plot(orb.x./host.Rvir,orb.vx./host.Vvir,'b.-',...
     'DisplayName','$X$');
 hold on 
 h(2)=plot(orb.y./host.Rvir,orb.vy./host.Vvir,'r.-',...
     'DisplayName','$Y$');
 grid
 xlabelmine('$x/R_\mathrm{vir}$')
 ylabelmine('$v/V_{\mathrm{vir}}$')
 hl=legend(h);
 set(hl,'Interpreter','latex','fontsize',14,'location','SouthWest');
 
 %% plot energy and AM conservation 
 figure
 h(1)=plot(tim,orb.energy(3,:)./orb.energy(3,1)-1,'b.-',...
     'DisplayName','$\delta E$');
 hold on 
 h(2)=plot(tim,orb.am./orb.am(1)-1,'r.-',...
     'DisplayName','$\delta L$');
 grid
 xlabelmine('$t\,[\mathrm{Gyr}]$')
 ylabelmine('$\delta E,\,\delta L$')
 hl=legend(h);
 set(hl,'Interpreter','latex','fontsize',14,'location','SouthWest');
 
 
 %% plot timestep 
 figure
plot(tim(2:end),diff(tim),'b.-');
 grid
 xlabelmine('$t\,[\mathrm{Gyr}]$')
 ylabelmine('$dt\,[\mathrm{Gyr}]$')
 hl=legend(h);
 set(hl,'Interpreter','latex','fontsize',14,'location','SouthWest');
 
 