function [ v1,d1] = FindEigpairs(A,K)
%FINDEIGPAIRS Summary of this function goes here
%   Detailed explanation goes here
[v,d]=eigs(A,K,'lr');

d=diag(d);

loc=(d==real(d)&d>0);
v=v(:,loc);
[~,L]=size(v);

v1=[];
d1=[];
thr=10^(-15);

for i=1:L
 v(v(:,i)<=thr&v(:,i)>=-thr,i)=0;

    
    if v(:,i)>=0
    v1=[v1 v(:,i)];
    d1=[d1 d(i)];
    end;

    if v(:,i)<=0
    v1=[v1 v(:,i)];
    d1=[d1 d(i)];
    end;
   
i
end;

if isempty(v1)
     warning('No strictly positive spectrums!');
end;
    
    
clear v d

end