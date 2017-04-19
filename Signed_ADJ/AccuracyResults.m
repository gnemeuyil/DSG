function [ACC,Qr,ANG,DB ] = AccuracyResults( n )
%ACCURACYRESULTS Summary of this function goes here
%   Detailed explanation goes here

k0 = [200,180,170,150,140];
k=k0*600;
%p=(triu(ones(5,5)*0.1,1)+triu(rand(5)*0.2,1)).*triu(ones(5,1)*(10./k),1);
%p=(triu(ones(5,5)*0.2,1)).*triu(ones(5,1)*(10./k),1);
%p=p+p';
ACC=[];
Qr=[];
ANG=[];
DB=[];
for i=1:n
    warning('Checking accuracy:')
    i
    tic;[Ap,A, partition] = KpartiteGenerator2(k);toc;
    tic;[~, ~, db0, Q0, ang0, acc0] = AdjCluster(Ap+A,10,2,0,partition);toc;
    tic;[~, ~, ~, ~, Q, acc] = NormalCutS(A+Ap,10,2,0,partition);toc;
    
        ACC(i,1)=acc0(4);
        ACC(i,2)=acc(4);
        ACC(i,3)=nnz(Ap);
        ACC(i,4)=nnz(A);
        Qr(i,1)=Q0(4);
        Qr(i,2)=Q(4);
        ANG(i,:)=ang0(4,:);
        DB(i,:)=db0(4,:);
  


end;



end

