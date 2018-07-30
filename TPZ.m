function X = TPZ(eval_f,eval_u,p,x_start,t_start,t_stop, timestep,N)
X(:,1) = x_start;
t_current = t_start;
for i=1:ceil((t_stop-t_start)/timestep)
   dt = min(timestep, (t_stop-t_current));
   u = feval(eval_u, t_current);
   f = feval(eval_f, X(:,i), u, p);
   X(:,i+1)=(eye(N)-0.5*dt*p.A)\(X(:,i)+0.5*dt*(f+p.B*u));
   t_current = t_current + dt;
%    plot(X);
%    pause(.2);
end
