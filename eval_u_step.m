function u = eval_u_step(t)
% generates the value of the input at time t
% corresponding to a step of magnitude a 
%
% u = eval_u_step(t);
 
if t <0
   u = 0;
else
   u = 10^6;
end
