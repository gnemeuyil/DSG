function [ L ] = UnNormLap( A )
%LAPLACIAN Summary of this function goes here
% Detailed explanation goes here
row = sum(A,2);
n = length(row);
D = spdiags(row(:),0,n,n);
L = D-A;
L = inv(L);
end


