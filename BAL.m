function [A_1, B_1, C_1, nodes_1] = BAL(A,B,C,D)
sys = ss(A, B, C, D); %convert to state-space mode
GRED = balancmr(sys); %when you run this it will pop up a graph and ask you the desired order. For heat bar just chose 1 or 2 as 2 singular values dominate the graph

%returns reduced model
A_1 = GRED.A; 
B_1 = GRED.B;
C_1 = GRED.C;

nodes_1 = length(C_1);
end