function [ L ] = laplacian( A )
%LAPLACIAN Summary of this function goes here
% Detailed explanation goes here
A=A|A';
A=double(A);
row = sum(A,2);
row = sqrt(row);
one = ones(size(row,1), 1);
% D^(-1/2)
row = sqrt(one./row);
n = length(row);
D = spdiags(row(:),0,n,n);
L =eye(n)-D*A*D;
end


