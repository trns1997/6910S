% %% Problem Setup
% clear all;
% close all;
% clc;

addpath(genpath([pwd,'/cbrewer']));
CT=cbrewer('qual','Set1',9);

N=50; % Nodes
[A,B,C,E,D]=virus(N); % Setup viral model
% [A,B,C,E,D] = HB(N);
p.A   = A; % 1x1 struct containing A and B
p.B   = B;
x_start = zeros(N,1); % initial conditions
% x_start(1,1) = 10^3; % initial conditions
t_start = 0; % start time
t_stop  = 1; % stop time
timestep = 1e-4; % time step 
dt1=0.0001; % variable time step for BE 1st part
dt2=0.1; % variable time step for BE 2nd part

%test input evaluation functions
% u  = eval_u_step(t_start);
eval_u = 'eval_u_step';
u = feval(eval_u,t_start); % evaluate u at t_start
%u = 0;

%eval_u = 'eval_u_impulse';
%u = feval(eval_u,t_input);

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

[X] = FE(eval_f,eval_u,p,x_start,t_start,t_stop,timestep); % returns a matrix of x evaluated using FE for all time points

[Y] = C*X; % evaluates y for all time points

%% Viz
n=5;
cntrsx=(1:100:100*n)';
cntrsy=14*ones(n,1);
strns=[1 12 25 37 50];

%%

i=2;
factor=4.5;
fig_prop(10,8);
for j=1:5
    hold on
    k=strns(1,j);
    pos=[cntrsx(j,1) cntrsy(j,1)-(X(k,i))/factor 2*(X(k,i))/factor 2*(X(k,i))/factor];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',CT(j,:));
    axis equal;
    xlim([-10 450]);
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
axis off
    text(cntrsx(j,1)+(X(k,i))/factor-5,cntrsy(j,1)+2*(X(k,i))/factor,num2str(j),'FontSize',20);
end

for j=1:4
    k=strns(1,j);
    l=k+1;
    fill([cntrsx(j,1)+2*(X(k,i))/factor,cntrsx(j,1)+2*(X(k,i))/factor,cntrsx(j+1,1),cntrsx(j+1,1)],[cntrsy(j,1)-A(k,l)/80,cntrsy(j,1)+A(k,l)/80,cntrsy(j,1)+A(l,k)/80,cntrsy(j,1)-A(l,k)/80],CT(6,:));
    axis equal;
end

pause(0.1)

%%

n=5;
cntrsx=(1:2000:2000*n)';
cntrsy=14*ones(n,1);
cntrs=[cntrsx,cntrsy];
strns=[1 12 25 37 50];

%%

factor=4.5;
factor1=4;
fig_prop(10,8);
for i=3:10:100
    clf;
    for j=1:5
        hold on
        k=strns(1,j);
        pos=[cntrsx(j,:) cntrsy(j,:)-(X(k,i))/factor 2*(X(k,i))/factor 2*(X(k,i))/factor];
        rectangle('Position',pos,'Curvature',[1 1],'FaceColor',CT(j,:));
        axis equal;
        xlim([-1000 14000]);
        ylim([-4000 4000]);
%         set(gca,'xtick',[]);
%         set(gca,'ytick',[]);
axis off
    end
    
    for j=1:4
        hold on
        k=strns(1,j);
        l=k+1;
        fill([cntrsx(j,:)+2*(X(k,i))/factor,cntrsx(j,:)+2*(X(k,i))/factor,cntrsx(j+1,:),cntrsx(j+1,:)],[cntrsy(j,:)-A(k,l)/factor1,cntrsy(j,:)+A(k,l)/factor1,cntrsy(j,:)+A(l,k)/factor1,cntrsy(j,:)-A(l,k)/factor1],CT(6,:));
        axis equal;
        xlim([-1000 14000]);
        ylim([-4000 4000]);
    end
    
    pause(0.5)
end

%% 
% 
% pts = 500;
% for t=1:pts
%     stem(X(:,t));
%     ylim([0,3.5*10^4]);
%     pause(0.01)
% end
% 
% 
% %% Eigenvalue Decomposition
% [A_evd,B_evd,C_evd,x1]= Modal_Analysis(A,B,C,x_start,N);
% %[A_evd, B_evd, C_evd, nodes_evd, V_evd] = EVD(A,B,C);
% x_start_evd = zeros(1,1);
% p_evd.A   = A_evd;
% p_evd.B   = B_evd;
% f_evd = eval_f_linear(x_start_evd,u,p_evd);
% eval_f_evd = 'eval_f_linear';
% 
% [X_evd] = FE(eval_f_evd,eval_u,p_evd,x_start_evd,t_start,t_stop, timestep);
% [Y_evd] = C_evd*X_evd;

%% Speed FE vs BE vs TPZ

fig_prop(10,8);
hold on;
b=bar(exec_time');
set(gca,'YScale','log');
xlim([0 4]);
for i=1:4
    b(i).FaceColor=CT(i,:);
end
xticks(1:5);
xticklabels({'10^{-4}','10^{-5}','10^{-6}'});
set(gca,'FontSize',15);
xlabel('\Deltat');
ylabel('Execution Time');
title('Comparison of Execution Time vs \Deltat for numerical integration methods');
legend('FE','BE','fBE','Tr','Location','SouthOutside','Orientation','Horizontal');
print('delta_t_exec_time','-dpng','-r300');

%% Error FE vs BE vs TPZ
% [x_ref] = TPZ(eval_f,eval_u,p,x_start,t_start,t_stop,10^-7,N);
% [y_ref] = C*x_ref;
% 
% error(1,1)= sum(abs(y_ref(1:1000:end) - y_fe_4))/sum(abs(y_ref(1:1000:end)));
% error(2,1)= sum(abs(y_ref(1:1000:end) - y_be_4))/sum(abs(y_ref(1:1000:end)));
% error(3,1)= sum(abs(y_ref(1:10000:end) - y_fbe_4(1:end-1)))/sum(abs(y_ref(1:1000:end)));
% error(4,1)= sum(abs(y_ref(1:1000:end) - y_tr_4))/sum(abs(y_ref(1:1000:end)));
% error(1,2)= sum(abs(y_ref(1:100:end) - y_fe_5(1:end-1)))/sum(abs(y_ref(1:100:end)));
% error(2,2)= sum(abs(y_ref(1:100:end) - y_be_5(1:end-1)))/sum(abs(y_ref(1:100:end)));
% error(3,2)= sum(abs(y_ref(1:1000:end) - y_fbe_5(1:end-1)))/sum(abs(y_ref(1:100:end)));
% error(4,2)= sum(abs(y_ref(1:100:end) - y_tr_5(1:end-1)))/sum(abs(y_ref(1:100:end)));
% error(1,3)= sum(abs(y_ref(1:10:end) - y_fe_6))/sum(abs(y_ref(1:10:end)));
% error(2,3)= sum(abs(y_ref(1:10:end) - y_be_6))/sum(abs(y_ref(1:10:end)));
% error(3,3)= sum(abs(y_ref(1:100:end) - y_fbe_6(1:end-1)))/sum(abs(y_ref(1:10:end)));
% error(4,3)= sum(abs(y_ref(1:10:end) - y_tr_6))/sum(abs(y_ref(1:10:end)));

fig_prop(10,8);
hold on;
b=bar(error');
set(gca,'YScale','log');
xlim([0 4]);
for i=1:4
    b(i).FaceColor=CT(i,:);
end
xticks(1:5);
xticklabels({'10^{-4}','10^{-5}','10^{-6}'});
set(gca,'FontSize',15);
xlabel('\Deltat');
ylabel('Error');
title('Comparison of Error vs \Deltat for numerical integration methods');
legend('FE','BE','fBE','Tr','Location','SouthOutside','Orientation','Horizontal');
print('error_exec_time','-dpng','-r300');


%% SVD
[X] = FE(eval_f,eval_u,p,x_start,t_start,t_stop,timestep);
[A_svd, B_svd, C_svd, nodes_svd, V_svd] = singVal(X,A,B,C,t_start,t_stop,timestep);
x_start_svd = V_svd'*x_start;
p_svd.A   = A_svd;
p_svd.B   = B_svd;
f_svd = eval_f_linear(x_start_svd,u,p_svd);
eval_f_svd = 'eval_f_linear';
[X_svd] = FE(eval_f_svd,eval_u,p_svd,x_start_svd,t_start,t_stop, timestep);
[Y_svd] = C_svd*X_svd;

%% Balance
[A_bal, B_bal, C_bal, nodes_bal] = BAL(A,B,C,D);
x_start_bal = zeros(nodes_bal,1);
p_bal.A   = A_bal;
p_bal.B   = B_bal;
f_bal = eval_f_linear(x_start_bal,u,p_bal);
eval_f_bal = 'eval_f_linear';
[X_bal] = FE(eval_f_bal,eval_u,p_bal,x_start_bal,t_start,t_stop, timestep);
[Y_bal] = C_bal*X_bal;

%% Krylov
[A_kryl, B_kryl, C_kryl, nodes_kryl, V_kryl] = KRYL(A,B,C);
x_start_kryl = V_kryl'*x_start;
p_kryl.A   = A_kryl;
p_kryl.B   = B_kryl;
f_kryl = eval_f_linear(x_start_kryl,u,p_kryl);
eval_f_kryl = 'eval_f_linear';
[X_kryl] = FE(eval_f_kryl,eval_u,p_kryl,x_start_kryl,t_start,t_stop, timestep);
[Y_kryl] = C_kryl*X_kryl;

%% Plot MOR comparison
fig_prop(10,6);
hold on;
plot(Y,'color',CT(1,:),'LineWidth',2);
plot(ey_svd_2,'color',CT(2,:),'LineWidth',2);
plot(ey_bal_2,'color',CT(3,:),'LineWidth',2);
plot(ey_kryl_2,'color',CT(4,:),'LineWidth',2);
xlim([t_start 1000]);

set(gca,'FontSize',15);
xlabel('time(s)');
ylabel('Output');
title('Error Performance of MOR methods(q = 2)');
% legend('Orig','SVD','KSA','Location','SouthOutside','Orientation','Horizontal');
legend('Orig','SVD','TBR','KSA','Location','SouthOutside','Orientation','Horizontal');
print('MOR_comp','-dpng','-r300');

%% zoom
fig_prop(10,6);
hold on;
plot(Y,'color',CT(1,:),'LineWidth',2);
plot(ey_svd_2,'o','color',CT(2,:),'LineWidth',2);
plot(ey_bal_2,'+','color',CT(3,:),'LineWidth',2);
plot(ey_kryl_2,'x','color',CT(4,:),'LineWidth',2);
xlim([465 515]);
ylim([1.173*10^4 1.218*10^4]);

set(gca,'FontSize',15);
xlabel('time(s)');
ylabel('Output');
title('Error Performance of MOR methods (q = 2)');
% legend('Orig','SVD','KSA','Location','SouthOutside','Orientation','Horizontal');
legend('Orig','SVD','TBR','KSA','Location','SouthOutside','Orientation','Horizontal');
print('MOR_comp','-dpng','-r300');


%% dT sweep
% exec_time=zeros(4,3); % 4 methods, 3 values of dT
% for i=4:6
%     dT = 10^(-i);
%     
%     [x_fe] = FE(eval_f,eval_u,p,x_start,t_start,t_stop,dT);
%     [y_fe] =  C*x_fe;
%     f = @() FE(eval_f,eval_u,p,x_start,t_start,t_stop,dT);
%     exec_time(1,i-3)=timeit(f);
%     
%     [x_be] = BE(eval_f,eval_u,p,x_start,t_start,t_stop,dT,N);
%     [y_be] = C*x_be;
%     f = @() BE(eval_f,eval_u,p,x_start,t_start,t_stop,dT,N);
%     exec_time(2,i-3)=timeit(f);
%     
%     [x_fbe] = BE_fast(eval_f,eval_u,p,x_start,t_start,t_stop,dT,N,dT,0.1);
%     [y_fbe] = C*x_fbe;
%     f = @() BE_fast(eval_f,eval_u,p,x_start,t_start,t_stop,dT,N,dT,0.1);
%     exec_time(3,i-3)=timeit(f);
%     
%     [x_tr] = TPZ(eval_f,eval_u,p,x_start,t_start,t_stop,dT,N);
%     [y_tr] = C*x_tr;
%     f = @() TPZ(eval_f,eval_u,p,x_start,t_start,t_stop,dT,N);
%     exec_time(4,i-3)=timeit(f);
%     
%     eval(['y_fe_' num2str(i) '= y_fe;']);
%     eval(['y_be_' num2str(i) '= y_be;']);
%     eval(['y_fbe_' num2str(i) '= y_fbe;']);
%     eval(['y_tr_' num2str(i) '= y_tr;']);
%     
%     eval(['x_fe_' num2str(i) '= x_fe;']);
%     eval(['x_be_' num2str(i) '= x_be;']);
%     eval(['x_fbe_' num2str(i) '= x_fbe;']);
%     eval(['x_tr_' num2str(i) '= x_tr;']);
%     
% end

%% Speed MOR N --> q
% s_orign = @() FE(eval_f,eval_u,p,x_start,t_start,t_stop,timestep);
% exec_time(1,1)=timeit(s_orign);
% 
% %500
% s_svd = @() FE(eval_f_svd,eval_u,p_svd,x_start_svd,t_start,t_stop, timestep);
% exec_time(1,2)=timeit(s_svd);
% s_bal = @() FE(eval_f_bal,eval_u,p_bal,x_start_bal,t_start,t_stop, timestep);
% exec_time(1,3)=timeit(s_bal);
% s_kryl = @() FE(eval_f_kryl,eval_u,p_kryl,x_start_kryl,t_start,t_stop, timestep);
% exec_time(1,4)=timeit(s_kryl);
% 
% %250
% s_svd = @() FE(eval_f_svd,eval_u,p_svd,x_start_svd,t_start,t_stop, timestep);
% exec_time(1,5)=timeit(s_svd);
% s_bal = @() FE(eval_f_bal,eval_u,p_bal,x_start_bal,t_start,t_stop, timestep);
% exec_time(1,6)=timeit(s_bal);
% s_kryl = @() FE(eval_f_kryl,eval_u,p_kryl,x_start_kryl,t_start,t_stop, timestep);
% exec_time(1,7)=timeit(s_kryl);
% 
% %47
% s_svd = @() FE(eval_f_svd,eval_u,p_svd,x_start_svd,t_start,t_stop, timestep);
% exec_time(1,8)=timeit(s_svd);
% s_bal = @() FE(eval_f_bal,eval_u,p_bal,x_start_bal,t_start,t_stop, timestep);
% exec_time(1,9)=timeit(s_bal);
% s_kryl = @() FE(eval_f_kryl,eval_u,p_kryl,x_start_kryl,t_start,t_stop, timestep);
% exec_time(1,10)=timeit(s_kryl);
% 
% %2
% s_svd = @() FE(eval_f_svd,eval_u,p_svd,x_start_svd,t_start,t_stop, timestep);
% exec_time(1,11)=timeit(s_svd);
% s_bal = @() FE(eval_f_bal,eval_u,p_bal,x_start_bal,t_start,t_stop, timestep);
% exec_time(1,12)=timeit(s_bal);
% s_kryl = @() FE(eval_f_kryl,eval_u,p_kryl,x_start_kryl,t_start,t_stop, timestep);
% exec_time(1,13)=timeit(s_kryl);

% fig_prop(10,8);
% b=bar(exec(2:end,:));
% for i=1:3
%     b(i).FaceColor=CT(i,:);
% end
% xticklabels({'500','250','47','2'});
% % set(gca,'YScale','log');
% grid on;
% set(gca,'FontSize',15);
% xlabel('Reduced Model Order (q)');
% ylabel('Time(s)');
% title('Comparison of Execution Time vs Model Order for MOR methods');
% legend('SVD','TBR','KSA','Location','SouthOutside','Orientation','Horizontal');
% text(3,1,['Original: ',num2str(exec_time_mor(1,1)),'(s)'],'FontSize',20);
% print('speed_MOR','-dpng','-r300');

fig_prop(10,8);
b=bar([exec_time_mor(1,1);exec(2:end,1)]);
b(1).FaceColor=CT(2,:);
xticklabels({'Orig','500','250','47','2'});
set(gca,'YScale','log');
ylim([3e-1 1e1]);
grid on;
set(gca,'FontSize',15);
xlabel('Reduced Model Order (q)');
ylabel('Time(s)');
title('Comparison of Execution Time vs Model Order');
print('speed_MOR','-dpng','-r300');

%% Error MOR N --> q
% ey_origin = Y;
% 
% %10
% ey_svd_10 = Y_svd;
% error(1,1) = sum(abs(ey_origin - ey_svd_500))/sum(abs(ey_origin));
% ey_bal_10 = Y_bal;
% error(1,2) = sum(abs(ey_origin - ey_bal_500))/sum(abs(ey_origin));
% ey_kryl_10 = Y_kryl;
% error(1,3) = sum(abs(ey_origin - ey_kryl_500))/sum(abs(ey_origin));
% 
% %5
% ey_svd_5 = Y_svd;
% error(2,1) = sum(abs(ey_origin - ey_svd_250))/sum(abs(ey_origin));
% ey_bal_5 = Y_bal;
% error(2,2) = sum(abs(ey_origin - ey_bal_250))/sum(abs(ey_origin));
% ey_kryl_5 = Y_kryl;
% error(2,3) = sum(abs(ey_origin - ey_kryl_250))/sum(abs(ey_origin));
% 
% %2
% ey_svd_2 = Y_svd;
% error(3,1) = sum(abs(ey_origin - ey_svd_2))/sum(abs(ey_origin));
% ey_bal_2= Y_bal;
% error(3,2) = sum(abs(ey_origin - ey_bal_2))/sum(abs(ey_origin));
% ey_kryl_2 = Y_kryl;
% error(3,3) = sum(abs(ey_origin - ey_svd_2))/sum(abs(ey_origin));

fig_prop(10,8);
b = bar(error_mor);
for i=1:3
    b(i).FaceColor=CT(i,:);
end
hold on
xticklabels({'10','5','2'});
grid on;
set(gca,'YScale','log');
set(gca,'FontSize',15);
xlabel('Reduced Model Order (q)');
ylabel('Error');
title('Comparison of Error vs Model Order for MOR methods');
legend('SVD','TBR','KSA','Location','SouthOutside','Orientation','Horizontal');
print('error_MOR','-dpng','-r300');

%% 
% 
% fig_prop(10,8);
% plot(y_fe_4,'LineWidth',2,'color',CT(1,:));
% hold on;
% set(gca,'FontSize',15);
% dummy=1.216e4*ones(1,1e3);
% plot(dummy,'k-.','LineWidth',2);
% xlim([0 1e3]);
% line([168 168],[0 14e3],'LineWidth',2,'LineStyle','--');
% text(40,13000,'Transient','FontSize',15);
% text(500,13000,'SS','FontSize',15);
% text(100,6000,'\Deltat_1','FontSize',15);
% text(500,6000,'\Deltat_2','FontSize',15);
% grid on;
% xlabel('Time (s)');
% ylabel('y(t)');
% title('System response for choosing \Deltat for fBE');
% print('fBE','-dpng','-r300');