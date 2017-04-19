function B = PLUS2(alpha,dmax,Nd)
%alpha is usually 1~1.5
%dmax is the max degree
%Nd is the number of bins

logdi = (0:Nd) * log(dmax) / Nd;
di = unique(round(exp(logdi)));
logni = alpha * (log(dmax) - log(di));
ni = round(exp(logni));
v = PowerLawEdges(di,ni);
%v(randperm(length(v)))
%idx=speye(size(max(v)));
idx=randperm(max(v));
B=triu(sparse(v,v(randperm(length(v))),1,max(v),max(v)),1);
B=B+B';
B=B(idx,idx);
