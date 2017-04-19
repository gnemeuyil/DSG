function [HB,D2] = WNormLap( A )
%WNORMLAP Summary of this function goes here
%   Detailed explanation goes here
[n,~]=size(A);

row = sum(A,2);
row2 = sqrt(row);
one = ones(n,1);
row3 = one./row2;

D = spdiags(row(:),0,n,n);
D2= spdiags(row3(:),0,n,n);
%HB = 0.5*D2*(2*D-A-A')*D2;
HB = spdiags(one(:),0,n,n)-0.5*D2*(A+A')*D2;



end

