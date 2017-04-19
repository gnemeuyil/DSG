function R = RandUS_Di_Neg(n,k,p)
%   Modified random asymmetrix matrix of size n with density k
%   The result will be a square matrix with Perron Frobenius property by 
%   making nodes strongly connected. The entries directely above the trace
%   are set to 1 to make the associated graph a SCC.
%   Negative edges are added such that -/+ is of percentage p

%generate random sequence of indices for upper triagular matrix
%I=randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]);

%convert to x-y coordinates accordin to the random indeices generated

%Upper half
[ xu,yu ] = TriuIdx2Coord( randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]));
%lower half
[ xl,yl ] = TriuIdx2Coord( randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]));

R=sparse([xu; yl],[yu; xl],1,n,n);
R=double(logical(R));
[x,y]=find(R);


I=sub2ind([n, n],x,y);% convert edge subscript to index
Eye_ind=sub2ind([n, n],1:n,1:n);% find the index of the diagonal elements
IND=1:n^2;% this is the index of all the possible edge entries
Edges = setdiff(IND,[I', Eye_ind]);% this is available edges indices without repeating the sub graph and diagonal elements

S = datasample(Edges,round(p*length(x)),'Replace',false);%sample some without replacement to conduct the test
R(S)=-1;


%set triu(1) to 1 to maintain SCC
R(n,1)=1;
i=1:(n-1);
R(i*(n+1))=1;


end

