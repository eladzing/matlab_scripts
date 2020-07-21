function res = rk4_orbitIntegration(IC,dt,tmax,dtmin,rhs,hostHalo,varargin)

%RK4_ORBITINTEGRATION - use 4th order runge-kutta to integrate an orbit
% hostHalo is a object of a halo class (eg. NFW).

%% setup
epsTot=1e-3;
epsStep=1e-4;
epsLow=1e-6;

cntMax=1e7;


%condFlag=false;
% rCondOper='';
% rCondVal=[];
cond='';
%% read conditions 

i=1;
while i<=length(varargin)
    switch(lower(varargin{i}))
%         case {'rcond','condr'}
%             i=i+1;
%             rCondOper=varargin{i};
%             i=i+1;
%             rCondVal=varargin{i};
%             
%             rCondFlag=true;
        case {'condition'  'cond'}
            i=i+1;
            cond=varargin{i};
            
        otherwise
            error('rk4_orbitIntegration - illegal argument: %s',varargin{i})
    end
    i=i+1;
end




%% set initial conditions


len=min(100.*ceil(tmax/dt),cntMax);

x=zeros(1,len);
y=x;
vx=x;
vy=x;
t=x;
nrg=zeros(3,len);
rad=x;
am=x;

x(1)=IC.x;
y(1)=IC.y;
vx(1)=IC.vx;
vy(1)=IC.vy;
t(1)=0;
rad(1)=sqrt(x(1)^2+y(1)^2);
nrg(1,1)=0.5*(vx(1)^2+vy(1)^2);
nrg(2,1)=potential(hostHalo,rad(1),'kpc');
nrg(3,1)=nrg(1,1)+nrg(2,1);
am(1)= x(1)*vy(1)-y(1)*vx(1);
%% begin integration
cnt=1;
%hb = waitbar(0,'Processing orbit...');
%hf=figure;
while t(cnt)<tmax && cnt<cntMax
    cnt=cnt+1;
    % preform 1 step
    
%     if mod(cnt,100)==0
%         waitbar(t(cnt-1)./tmax,hb,sprintf('dt=%g',dt))
%     end
    
    [x(cnt),y(cnt),vx(cnt),vy(cnt)] = orbits.rk4_step(x(cnt-1),y(cnt-1),vx(cnt-1),vy(cnt-1),dt,rhs,hostHalo);
    t(cnt)=t(cnt-1)+dt;
    
    %% check error & adapt timestep
    
    % calculate energy
    rad(cnt)=sqrt(x(cnt)^2+y(cnt).^2);
    nrg(1,cnt)=0.5*(vx(cnt)^2+vy(cnt)^2);
    nrg(2,cnt)=potential(hostHalo,rad(cnt),'kpc');
    am(cnt)= x(cnt)*vy(cnt)-y(cnt)*vx(cnt);
    
    nrg(3,cnt)=nrg(1,cnt)+nrg(2,cnt);
    % check error and adap timestep
    dETot=abs(nrg(3,cnt)/nrg(3,1)-1);
    dLTot=abs(am(cnt)/am(1)-1);
    dE=abs(nrg(3,cnt)/nrg(3,cnt-1)-1);
    dL=abs(am(cnt)/am(cnt-1)-1);

    
    %% check ext conditions 
    
    if ~isempty(cond)
       %cond=['rad(cnt) ' rCondOper ' ' num2str(rCondVal)];
       
       if eval(cond)
           break
       end
    end
    
%     if dETot>epsTot || dLTot>epsTot  || dE>epsStep || dL>epsStep  % if energy or AM is not conserved reduce timestep and rerun step
%         
%         dt=0.5*dt;
%         cnt=cnt-1;
%         
%         if dt<dtmin % don't go below minimal timestep
%             dt=dtmin;
%             cnt=cnt+1;
%         end
%         
%     elseif dE<epsLow && dL<epsLow
%         dt=2*dt;
%     end
  
end
% close(hb) 

res.x=x(1:cnt);
res.y=y(1:cnt);
res.rad=rad(1:cnt);
res.vx=vx(1:cnt);
res.vy=vy(1:cnt);
res.t=t(1:cnt);
res.energy=nrg(:,1:cnt);
res.am=am(1:cnt);
res.vel=sqrt(vx(1:cnt).^2+vy(1:cnt).^2);




end

