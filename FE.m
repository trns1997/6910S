function X = FE(eval_f,eval_u,p,x_start,t_start,t_stop, timestep)
% uses Forward Euler to simulate states model dx/dt=f(x,u,p)
% from state x_start at time t_start
% until time t_stop, with time intervals timestep
% eval_f is a string including the name of the function that evaluates f(x,u,p)
% eval_u os a string including the name of the funciton that evaluates u(t)
% 
% X = FE(eval_f,p,x_start,t_start,t_stop, dt)

X(:,1) = x_start;
t_current = t_start;
for i=1:ceil((t_stop-t_start)/timestep)
   dt = min(timestep, (t_stop-t_current));
   u = feval(eval_u, t_current);
   f = feval(eval_f, X(:,i), u, p);
   X(:,i+1)= X(:,i) +  dt * f;
%    X(:,i+1)=X(:,i+1)./sum(X(:,i+1));
   t_current = t_current + dt;
end

% plot(X);
end