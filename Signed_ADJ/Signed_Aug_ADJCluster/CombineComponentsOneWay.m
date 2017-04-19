function [G ] = CombineComponentsOneWay( A,B,p,IO)
%COMBINECOMPONENTSONEWAY Summary of this function goes here
% 
% IO==1, edges are all from A to B
% IO==2, edges are all from B to A
sa=size(A,1);
sb=size(B,1);

G=zeros(sa+sb, sa+sb);
G(1:sa,1:sa)=A;
G(sa+1:sa+sb,sa+1:sa+sb)=B;

NA=ones(sa,sb);
NB=ones(sb,sa);
PA=NA.*p;
PB=NB.*p;
A2B=binornd(NA,PA);
B2A=binornd(NB,PB);

if IO==1
   G(1:sa,sa+1:sa+sb)=A2B; 
end

if IO==2
   G(sa+1:sa+sb,1:sa)=B2A;
end


end

