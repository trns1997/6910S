function X = BE_fast(eval_f,eval_u,p,x_start,t_start,t_stop, timestep,N,dt1,dt2)
% uses Backward Euler to simulate states model dx/dt=f(x,u,p)
% from state x_start at time t_start
% until time t_stop, with time intervals timestep
% eval_f is a string including the name of the function that evaluates f(x,u,p)
% eval_u os a string including the name of the funciton that evaluates u(t)
% 
% X = FE(eval_f,p,x_start,t_start,t_stop, dt)
t1 = length(0:dt1:0.1); %get indices based on dT1 which ranges from 10^(i = 4:6)
t2 = length(0.1:dt2:(t_stop-0.1)); %get indices based on dT2 which is 0.1 in the main when you call this function

X(:,1) = x_start;
t_current = t_start;
for i=1:t1
   dt = min(timestep, (t_stop-t_current));
   u = feval(eval_u, t_current);
   f = feval(eval_f, X(:,i), u, p);
   X(:,i+1)=(eye(N)-dt1*p.A)\(X(:,i)+dt1*p.B*u);
   t_current = t_current + dt1;
%    plot(X);
%    pause(.2);
end
for i=t1:(t1+t2)
   dt = min(timestep, (t_stop-t_current));
   u = feval(eval_u, t_current);
   f = feval(eval_f, X(:,i), u, p);
   X(:,i+1)=(eye(N)-dt2*p.A)\(X(:,i)+dt2*p.B*u);
   t_current = t_current + dt2;
%    plot(X);
%    pause(.2);
end
end
