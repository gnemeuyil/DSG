function R = RandUS(n,k)
%   Modified random symmetrix matrix of size n with density k
%   The result will be a sparse matrix

%generate random sequence of indices for upper triagular matrix
%I=randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]);

%convert to x-y coordinates accordin to the random indeices generated
[ x,y ] = TriuIdx2Coord( randi(n*(n-1)/2,[round(k*n*(n-1)/2),1]));
R=sparse(x,y,1,n,n);
R = R+R';







