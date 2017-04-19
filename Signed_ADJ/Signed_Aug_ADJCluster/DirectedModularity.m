function [ A ] = DirectedModularity( M, J )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%  J c
P=idx2lgc(J);

KIn = sum(M,1); % Sum of In weight with node i
KOut = sum(M,2); % Sum of Out weight with node i
SumTotIn = sum(M,1);
SumTotOut = sum(M,2);
SumIn = diag(M);
m = sum(sum(M));

delt=spase(size(J), size(J));
for i=1:max(J{2})
    C=J{2}==i;
    delt[C,:]=



end

