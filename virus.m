function [A,B,C,E,D] = virus(nodes)

E = eye(nodes);

A = 2500*(-eye(nodes)+diag(sqrt(2:nodes),1)/nodes+diag(2:nodes,-1)/nodes);

B = ones(nodes,1);
%B(nodes/2,1)=1;

% C = eye(nodes);

C = zeros(1, nodes);
C(nodes)=1;

D = 0;

end