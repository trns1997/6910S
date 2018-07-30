function [A_1, B_1, C_1, nodes_1, U1] = singVal(X,A,B,C,t_start,t_stop,timestep)
x_comp = X(1:end, 1:10:ceil((t_stop-t_start)/timestep)); %sample training points

[U,S] = svd(x_comp);

stem(diag(S(1:10,1:10))); %to determine which eigenvalues to keep
xlim([0 11]);
set(gca,'YScale','log');
xlabel('Singular values');
ylabel('Magnitude');
title('Plot of largest 10 singular values for MOR');
print('sing_val','-dpng','-r300');

U1 = U(1:end, 1:2); %keeps only the first column

B_h = U1'*B;
C_h = C*U1;
A_h = U1'*A*U1;

A_1 = A_h;
B_1 = B_h;
C_1 = C_h;
nodes_1 = length(C_1);
end
