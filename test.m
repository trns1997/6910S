%% Problem Setup
clear all;
close all;
clc;

N=50; % Nodes
[A,B,C,E,D]=HB(N); % Setup heat bar
p.A   = A; % 1x1 struct containing A and B
p.B   = B;
x_start = zeros(N,1); % initial conditions
t_start = 0; % start time
t_stop  = 1; % stop time
timestep = 1e-4; % time step 
dt1=0.0001; % variable time step for BE 1st part
dt2=0.1; % variable time step for BE 2nd part

%test input evaluation functions
% u  = eval_u_step(t_start);
eval_u = 'eval_u_step';
u = feval(eval_u,t_input); % evaluate u at t_start
%u = 0;

% test vector field evaluation functions
f = eval_f_linear(x_start,u,p); 
eval_f = 'eval_f_linear'; % evaluate f at t_start

%f = eval_f_SquareDiagonal(x_start,u,p)
%eval_f = 'eval_f_SquareDiagonal';
% f = feval(eval_f,x_start,u,p);

% test Jacobian evaluation functions
%Jf = eval_Jf_linear(x_start,u,p)
%eval_Jf = 'eval_Jf_linear';
%Jf = eval_Jf_SquareDiagonal(x_start,u,p)
%eval_Jf = 'eval_Jf_SquareDiagonal';
%Jf = feval(eval_Jf,x_start,u,p)

%% test FE function
[X] = FE(eval_f,eval_u,p,x_start,t_start,t_stop, timestep); % returns a matrix of x evaluated using FE for all time points
[Y] = C*X; % evaluates y for all time points
plot(Y); % plots y
hold on

%% Eigenvalue Decomposition
[A_evd, B_evd, C_evd, nodes_evd] = EVD(A,B,C);
x_start_evd = zeros(nodes_evd,1);
p.A   = A_evd;
p.B   = B_evd;
f_evd = eval_f_linear(x_start_evd,u,p);
eval_f_evd = 'eval_f_linear';

[X_evd] = FE(eval_f_evd,eval_u,p,x_start_evd,t_start,t_stop, timestep);
[Y_evd] = C_evd*X_evd;
plot(Y_evd);

%% SVD
[X] = FE(eval_f,eval_u,p,x_start,t_start,t_stop, timestep);
[A_svd, B_svd, C_svd, nodes_svd] = SVD(X,A,B,C,t_start,t_stop,timestep);
x_start_svd = zeros(nodes_svd,1);
p.A   = A_svd;
p.B   = B_svd;
f_svd = eval_f_linear(x_start_svd,u,p);
eval_f_svd = 'eval_f_linear';
[X_svd] = FE(eval_f_svd,eval_u,p,x_start_svd,t_start,t_stop, timestep);
[Y_svd] = C_svd*X_svd;
plot(Y_svd);

%% Balance
[A_bal, B_bal, C_bal, nodes_bal] = BAL(A,B,C,D);
x_start_bal = zeros(nodes_bal,1);
p.A   = A_bal;
p.B   = B_bal;
f_bal = eval_f_linear(x_start_bal,u,p);
eval_f_bal = 'eval_f_linear';
[X_bal] = FE(eval_f_bal,eval_u,p,x_start_bal,t_start,t_stop, timestep);
[Y_bal] = C_bal*X_bal;
plot(Y_bal);