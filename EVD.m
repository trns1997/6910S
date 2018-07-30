function [A_1, B_1, C_1, nodes_1, V] = EVD(A,B,C)

[V,D] = eig(A);

B_h = V\B;
C_h = C*V;
D_h = diag(D);

out_h = (B_h.*C_h')./D_h; % divided by the eigenvalue to clearly see the dominant pole. big out_h means eigenvalue very small 
stem(abs(out_h))

ind = abs(out_h) > 1.684*10^9; %to get the index of the smallest eigenvalue. Value determined from above plot.

C_h = C_h(ind);

B_h = B_h(ind);

D_h1 = D_h(ind);

A_1 = diag(D_h1);
B_1 = B_h;
C_1 = C_h;
nodes_1 = length(C_1);

end