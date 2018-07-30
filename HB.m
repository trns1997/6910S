function [A,B,C,E,D] = HB(nodes)

E = eye(nodes);

A = 2*eye(nodes)+diag(-1*ones(1,nodes-1),1)+diag(-1*ones(1,nodes-1),-1);
A(1,1) = 1;
A(nodes,nodes) = 1;
A = -((nodes^2))*A;

B = zeros(nodes,1);
B(20)=1;

C = zeros(1, nodes);
C(nodes)=1;

D = 0;

end
