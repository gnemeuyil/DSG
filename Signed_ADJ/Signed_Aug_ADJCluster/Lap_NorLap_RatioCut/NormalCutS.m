function [Jf, Cf, kf, db, Q, acc] = NormalCutS(A,K,L,c,Partition)

m = size(A,1);
if ~exist('K','var') && ~exist('L','var')&& ~exist('c','var')
    K = min(100, ceil(m/10));
    L = 2;
    c = 0;
end;

%d = sparse(diag(sum(abs(A))+1));
As=A+A';
d = sparse(1:m,1:m,sum(abs(As)).^(-0.5),m,m);
I=sparse(1:m,1:m,1,m,m);
N =I- d*As*d;
%d = sparse(diag(sum(abs(A))));
%N = d-A;
[V,D]=eigs(N,K,'SM');

db = zeros(K-1,3);
kf = zeros(3,1);
Q = zeros(K-1,1);
acc = zeros(K-1,1);
temp = [10000 0];
Cf = cell(3,1);
Df = cell(3,1);


for i = 2:K
    [Jt, Ct, dbt] = AdjCl(V(:,1:i),L,c);
    db(i-1,:) = dbt;
    if db(i-1,1)<temp(1)
        temp(1) = db(i-1,1);
        Jf{1} = Jt;
        Cf{1} = Ct;
        kf(1) = i;
    elseif temp(1) == db(i-1,1)
        %WARNING('Multiple Solution!)');
    end;
    
    P = idx2lgc(Jt);
    %Q(i-1) = Qfunction(A,P);
    Q(i-1) = SignQfunction(A,P);
    if Q(i-1)>temp(2)
        temp(2) = Q(i-1);
        Jf{2} = Jt;
        Cf{2} = Ct;
        kf(2) = i;
    elseif temp(2) == Q(i-1)
        WARNING('Multiple Solution!)');
    end;
    
    if exist('Partition','var') && size(P,2)-size(Partition,2)>=0;
        PartitionE=zeros(size(P));
        PartitionE(:,1:size(Partition,2))=Partition;
        acc(i-1) = PartitionAccuracy(P, PartitionE);
    end
end;


% figure;plot(2:K,db(:,1),'-o','MarkerSize',6,'MarkerFaceColor','b');
% xlabel('k');
% ylabel('Davies-Buldin Index');
% hold on;
% plot([2,K],[1,1],'k--');
% hold off;
% figure;errorbar(2:K,db(:,1),db(:,2)-db(:,1),db(:,1)-db(:,3));
% xlabel('k');
% ylabel('Davies-Buldin Index');
% xlim([2,K]);
% hold on;
% plot([2,K],[1,1],'k--');
% hold off;
% 
% [qx,qy] = max(Q(logical(db(:,1)<1)));
% figure;plot(2:K,Q,'-o','MarkerSize',6,'MarkerFaceColor','b');
% xlabel('k');
% ylabel('Modularity');
% hold on;
% plot(qy+1,qx,'ks','MarkerSize',12)
% hold off;
% 
% d = diag(D);
% [dx,dy]=sort(abs(d),'descend');
% figure;plot(1:K,d(dy),'-o','MarkerSize',6,'MarkerFaceColor','b');
% xlabel('k');
% ylabel('Eigenvalues');
% hold on;
% plot([1,K],[0,0],'k--');
% hold off;

