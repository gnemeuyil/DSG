function R = PLUS(n,varargin)
%node number n
%alpha is usually 1~1.5
%2.1~2.5 2.3 common
m = size(varargin,2);


if m==1
    d = varargin{1};
    if size(d,1)==1
        d = d';
    end;
elseif m==2
    xmin = varargin{1};
    alpha = varargin{2};
    d = sparse(round(xmin*(rand(n,1)).^(-1/(alpha-1))));
    %d = sparse(round(xmin*(rand(n,1)).^(-alpha)));
end;
%P = d*d'/sum(d);% 30 minutes for size of 92987 with n=20k
%index=find(d>0);
%sd=sum(d);
%A=sparse(X,Y,rand<=d(X).*d(Y)/sd);
%R = logical(A<=P);
k=nnz(d)
% censor first then multiply to reduce the amout of elements
% d=d/sqrt(sum(d));
% r=sqrt(rand(n,1));
% d=d(r<=d);


warning('get edge list')
tic;
[x,y,s]=find(triu(d*d',1));
toc;
warning('get sparse matrix');
tic;
R=sparse(x,y,rand(size(s,1),1)<=(s/sum(d)),n,n);
toc;
R = R + R';