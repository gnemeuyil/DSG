function R = UniUS(n)

A = sprand(n);

d = unidrnd(n-4,n,1)+3;
P = d*d'/sum(d);
R = logical(A<=P);
R = triu(R,1);
R = R + R';