function [ G,a2b,b2a ] = CombineComponents( A,B,p,p2)
% Combines two components with edges of probanility p
% binomial with N=1 is equivalent to bernoli 
% A is bigger component B is smaller component

sa=size(A,1);
sb=size(B,1);
% penalize p with respect to the smaller community size
P=sb/sa;
p=p*P
p2=p2*P


G=zeros(sa+sb, sa+sb);
G(1:sa,1:sa)=A;
G(sa+1:sa+sb,sa+1:sa+sb)=B;

NA=ones(sa,sb);
NB=ones(sb,sa);
PA=NA.*p;
PB=NB.*p2;

A2B=binornd(NA,PA);
B2A=binornd(NB,PB);
a2b=sum(sum(A2B))
b2a=sum(sum(B2A))
G(1:sa,sa+1:sa+sb)=A2B;
G(sa+1:sa+sb,1:sa)=B2A;



end

