function [A1,b1,c1,x1]= Modal_Analysis(A0,b0,c0,x0,P)
[V,D]=eig(A0);

% %%%Q: complex V???
% % D all negative and we consider constant unit input, therefore
% % eigvalue with large abs are discarded
%
[~,idx]=sort(b0.*c0'./real(diag(D)),'descend');
idx=idx(1:P);
Vtmp=[];
i=1;
D_h = diag(real(D));


while i<=P
   if isreal(V(:,idx(i)))
        Vtmp=[Vtmp,V(:,idx(i))];
        i=i+1;
   else
       if i<=P-1
           Vrtmp=real(V(:,idx(i)));
            Vitmp=imag(V(:,idx(i)));
            Vtmp=[Vtmp,Vrtmp,Vitmp];
        end
        i=i+2;
    end
end







% [~,P]=size(Vtmp);
% 
% Vtmp(:,1)=Vtmp(:,1)/norm(Vtmp(:,1),2);
% for i=2:1:P
%     for j=1:1:i-1
%         Vtmp(:,i)=Vtmp(:,i)-Vtmp(:,i)'*Vtmp(:,j)*Vtmp(:,j);
%     end
%     Vtmp(:,i)=Vtmp(:,i)/norm(Vtmp(:,i),2);
% end

A1=Vtmp'*A0*Vtmp;
% A1=D(idx(1:P),idx(1:P));
b1=Vtmp'*b0;
c1=c0*Vtmp;
x1=Vtmp'*x0;

out_h = b0.*c0'./real(D_h); % divided by the eigenvalue to clearly see the dominant pole. big out_h means eigenvalue very small 
stem(abs(out_h))

ind = abs(out_h) > 0.01*10^-3; %to get the index of the smallest eigenvalue. Value determined from above plot.

c1 = c1(ind);

b1 = b1(ind);

D_h1 = D_h(ind);

A1 = diag(D_h1);
b1 = b1;
c1 = c1;




end

