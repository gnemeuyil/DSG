function Q=SignQfunction(A,x)

% Q=Qfunction(A,x,type)
% A: adjacency matrix
% x: partition matrix
% type: if type=2, output normalized value. Default: unnormalized value.

n = size(A,1);
d = full(sum(A,2));
m = sum(d);
if nargin<2
    x = 0;
end
if isscalar(x)
    if x<1
        x = 1;
    end
    Q = eigs(A-d*d'/m,x,'la');
    Q = n*sum(Q)/x/2/m;
else
    x=logical(x);
    CmNo=size(x,2);
    Q=0;
    for i=1:CmNo
        b=sum(sum(A(x(:,i),x(:,i))));
        c=sum(sum(A(x(:,i),:)));
        Q=Q+b/m-(c/m)*(c/m);
    end
end

% % Q=Qfunction(A,x,type)
% % A: adjacency matrix
% % x: partition matrix
% % type: if type=2, output normalized value. Default: unnormalized value.
% 
% if min(min(A))==0
%     Q = Qfunction(A,x);
% else
%     P = double(logical(A>0));
%     N = double(logical(A<0));
%     m = sum(sum(N));
%     q1 = Qfunction(P,x);
%     if nargin<2
%     x = 0;
%     end
%     if isscalar(x)
%         q2 = 0;
%     else
%         x=logical(x);
%         CmNo=size(x,2);
%         q2=0;
%         for i=1:CmNo
%             for j = (i+1):CmNo
%             b=2*nnz(N(x(:,i),x(:,j)));
%             c1=nnz(N(x(:,i),:));
%             c2=nnz(N(x(:,j),:));
%             q2=q2+b/m-(c1/m)*(c2/m);
%             end;
%         end
%         
%         
%     end
%     Q = q1+q2;
% end;
% 
