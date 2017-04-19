function [ v1,d1,v2,d2,I] = ScreenEigpairs(A,K)
%FINDEIGPAIRS Summary of this function goes here
%   Detailed explanation goes here
[v,d]=eigs(A,K,'lr');

d=diag(d);

loc=(d==real(d)&d>0);
v=v(:,loc);
[~,L]=size(v);

v2=[];
d2=[];
thr=10^(-8);
I=[];
for i=1:L
    tempv= v(:,i);
   
 tempv(tempv<=thr&tempv>=-thr)=0;

    
    if tempv>=0
    v2=[v2 v(:,i)];
    d2=[d2 d(i)];
    I=[I i];
    end;

    if tempv<=0
    v2=[v2 v(:,i)];
    d2=[d2 d(i)];
    I=[I i];
    
    end;
   
i
end;

if isempty(v2)
     warning('No strictly positive spectrums!');
end;
    v1=v;
    d1=d;
    
clear v d

end