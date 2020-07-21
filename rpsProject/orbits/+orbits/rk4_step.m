function [x,y,vx,vy] = rk4_step(x0,y0,vx0,vy0,dt,rhs,host)
% RK4_STEP perform one 4th order runge-kutta step 
%   given the position & velocity at the beginning of the time step dt
%   advance the particle based on the force function given in rhs .
%   rhs is a function handle


%% calculate components 


xt=x0;
yt=y0;
vxt=vx0;
vyt=vy0;

[xDot1,yDot1,vxDot1,vyDot1]=rhs(xt,yt,vxt,vyt,host);


xt=x0+0.5*dt*xDot1;
yt=y0+0.5*dt*yDot1;
vxt=vx0+0.5*dt*vxDot1;
vyt=vy0+0.5*dt*vyDot1;

[xDot2,yDot2,vxDot2,vyDot2]=rhs(xt,yt,vxt,vyt,host);


xt=x0+0.5*dt*xDot2;
yt=y0+0.5*dt*yDot2;
vxt=vx0+0*0.5*dt*vxDot2;
vyt=vy0+0.5*dt*vyDot2;

[xDot3,yDot3,vxDot3,vyDot3]=rhs(xt,yt,vxt,vyt,host);

xt=x0+dt*xDot3;
yt=y0+dt*yDot3;
vxt=vx0+dt*vxDot3;
vyt=vy0+dt*vyDot3;

[xDot4,yDot4,vxDot4,vyDot4]=rhs(xt,yt,vxt,vyt,host);



%% advance 
x=x0 + dt/6*(xDot1 + 2*xDot2 + 2*xDot3+ xDot4); 
y=y0 + dt/6*(yDot1 + 2*yDot2 + 2*yDot3+ yDot4); 

vx=vx0 + dt/6*(vxDot1 + 2*vxDot2 + 2*vxDot3+ vxDot4); 
vy=vy0 + dt/6*(vyDot1 + 2*vyDot2 + 2*vyDot3+ vyDot4); 



end

