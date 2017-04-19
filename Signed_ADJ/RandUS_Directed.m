function R = RandUS_Directed(n,k)
%   Modified random asymmetrix matrix of size n with density k
%   The result will be a square matrix with Perron Frobenius property by 
%   making nodes strongly connected. The entries directely above the trace
%   are set to 1 to make the associated graph a SCC.

%generate random sequence of indices for upper triagular matrix
%I=randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]);

%convert to x-y coordinates accordin to the random indeices generated

%Upper half
[ xu,yu ] = TriuIdx2Coord( randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]));
%lower half
[ xl,yl ] = TriuIdx2Coord( randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]));

R=sparse([xu; yl],[yu; xl],1,n,n);
%set triu(1) to 1 to create SCC
R(n,1)=1;
i=1:(n-1);
R(i*(n+1))=1;
R=double(logical(R));

end

