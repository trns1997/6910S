function [A_1, B_1, C_1, nodes_1, V] = KRYL(A,B,C)
[L,U] = lu(A);
V(:,1) = B/norm(B);
i = 1;
j = 1;
% for i = 1:4
    y = L\V(:,i);
    V(:,i+1) = U\y;
%     for j = 1:i
        V(:, i+1) = V(:, i+1) - (V(:, i+1)'*V(:, j))*V(:, j); 
%         V(:, i+1) = V(:, i+1) - (V(:, i+1)/norm(V(:, i+1))); 
%     end
    V(:, i+1) = V(:, i+1)/norm(V(:, i+1));
% end

% rank(V);

B_h = V'*B;
C_h = C*V;
A_h = V'*A*V;

A_1 = A_h;
B_1 = B_h;
C_1 = C_h;
nodes_1 = length(C_1);
end