function [ v1,d1,v2,d2,I] = ScreenEigpairsSignedC2M(A,K)
%FINDEIGPAIRS Summary of this function goes here
%   Select real eigen pairs by probing up to K eigen pairs of A
%   v1 d1 is the original eigenpair
%   v2 d2 is the candidate eigenpair
%   I is the index for candifate eigenpair
%   Dealing with signed graphs by mapping complex eigenvectors into one
%   real vector containing sign(real part)*modulus
[v,d]=eigs(A,K,'lm');
d=diag(d);

%remove nonconvergent NaN entries
d=d(~isnan(real(d)));
v=v(:,~isnan(real(d)));


% remove conjugates
[~,ia1,~] = unique(abs(d));
d=d(ia1);
v=v(:,ia1);

% Sort according to modulus of eigenvalues descend
[~,I]=sort(abs(d),'descend');
d1=d(I);
v1=v(:,I);

% Find real eigenpairs
loc=(d1==real(d1));
idx=find(loc);
d=d1(loc);
v=v1(:,loc);
[~,L]=size(v);
 
% convert complex eigenvectors to moduli
Nloc=~loc;
v1(:,Nloc)=sign(real(v1(:,Nloc))).*abs(v1(:,Nloc));

% Find eigenpairs with positive eigenvectors
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
    I=[I idx(i)];
    end;

    if tempv<=0
    v2=[v2 v(:,i)];
    d2=[d2 d(i)];
    I=[I idx(i)];
    end; 
    i
end;

if isempty(v2)
     warning('No strictly positive spectrums!');
end;
  


clear v d

end
