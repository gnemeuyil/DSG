function [ Idx ] = KNei(ind, A,K,n)
%KNEIB Summary of this function goes here
%   Detailed explanation goes here
% K is the expected number of neighbors
% n is the max path length to search
% ind is the vertex/node to search for
for i=1:n
    
    kneigh = kneighbors(A,ind,i);
    Idx=union(Idx, kneigh);
    S=size(Idx);
    if S>=K
        i=n+1;
    end
end



end

