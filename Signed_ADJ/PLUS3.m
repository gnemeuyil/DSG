function B = PLUS3(alpha,dmax,Nd,k1,k2)
% This function is used to generate inter cluster edges.
% alpha is usually 1~1.5
% dmax is the max degree which should be about 3~5% of min(k1,k2)
% Nd is the number of bins
%

logdi = (0:Nd) * log(dmax) / Nd;
di = unique(round(exp(logdi)));
logni = alpha * (log(dmax) - log(di));
ni = round(exp(logni));
v = PowerLawEdges(di,ni);
%v(randperm(length(v)))
%idx=speye(size(max(v)));
%idx=randperm(max(v));
%B=triu(sparse(v,v(randperm(length(v))),1,max(v),max(v)),1);
%B=B+B';
v=v(v<=min(k1,k2));
B=sparse(v,v(randperm(length(v))),1,k1,k2);

B=B(randperm(max(k1)),randperm(max(k2)));
%B=B(idx,idx);


