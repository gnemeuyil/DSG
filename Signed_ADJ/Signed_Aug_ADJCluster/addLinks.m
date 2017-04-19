function [ G ,A2B, B2A] = addLinks( G,s,x1,x2,dx1,dx2,p1,p2 )
%ADDLINKS Summary of this function goes here
%   Detailed explanation goes here
sa=size(G,1)-s;

% penalize p with respect to the smaller community size
P=s/sa;
p1=p1*P
p2=p2*P


NA=ones(dx1,dx2);
NB=ones(dx2,dx1);
PA=NA.*p1;
PB=NB.*p2;

A2B=binornd(NA,PA);
B2A=binornd(NB,PB);
a2b=sum(sum(A2B))
b2a=sum(sum(B2A))
G(x1:x1+dx1-1,x2:x2+dx2-1)=A2B;
G(x2:x2+dx2-1,x1:x1+dx1-1)=B2A;


end

