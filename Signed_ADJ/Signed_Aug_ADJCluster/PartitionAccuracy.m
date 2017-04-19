function [tolM,M] = PartitionAccuracy(P, StdP)
% P: the partition matrix to be compared.
% StdP: the standard partition matrix as the benchmark.

M = double(P')*double(StdP);
Match = Hungarian(-M);
M = M(logical(Match));
tolM = sum(M)/nnz(StdP);
M = M./sum(StdP)';
