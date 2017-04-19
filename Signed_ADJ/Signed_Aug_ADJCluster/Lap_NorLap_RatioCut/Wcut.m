function [ Jf, Cf, kf, db, Q, acc, V, D] = Wcut( B,K,L,c,C,Partition)
%WCUT weighted cut that performs spectral projection using Hermitian of A
% associated with the degree of nodes

[HB,D2] = WNormLap( B );
A=HB;

m = size(A,1);
if ~exist('K','var') && ~exist('L','var')&& ~exist('c','var')
    K = min(100, ceil(m/10));
    L = 2;
    c = 0;
end;

%d = sparse(diag(sum(abs(A))+1));
%N = d^(-1)*A;
%d = sparse(diag(sum(abs(A))));
%N = d-A;
[V,DA]=eigs(A,K,'SA');

db = zeros(K-1,3);
kf = zeros(3,1);
Q = zeros(K-1,1);
acc = zeros(K-1,1);
temp = [10000 0];
Cf = cell(3,1);
Jf = cell(3,1);


for i = 2:K
    %[Jt, Ct, dbt] = AdjCl(D2*V(:,1:i),L,c);
    [Jt, dbt, ~, Ct] = WAdjCl(D2*V(:,1:i),L,c);
    db(i-1,:) = dbt;
    if db(i-1,1)<temp(1)
        temp(1) = db(i-1,1);
        Jf{1} = Jt;
        Cf{1} = Ct;
        kf(1) = i;
    elseif temp(1) == db(i-1,1)
        %WARNING('Multiple Solution!');
    end;
    
    P = idx2lgc(Jt);
    %Q(i-1) = Qfunction(A,P);
    Q(i-1) = SignQfunction(B,P);
    if Q(i-1)>temp(2)
        temp(2) = Q(i-1);
        Jf{2} = Jt;
        Cf{2} = Ct;
        kf(2) = i;
    elseif temp(2) == Q(i-1)
       temp(2)
    end;
    
    if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
        acc(i-1) = PartitionAccuracy(P, Partition);
    end
end;

V=V(:,1:kf(2));
D=DA(:,1:kf(2));
figure;
plot(2:K,(db(:,1)-mean(db(:,1)))/(max(db(:,1))-min(db(:,1))),'-o','MarkerSize',6,'MarkerFaceColor','b');
xlabel('k');
ylabel('Davies-Buldin Index');
hold on;
plot([2,K],[1,1],'k--');
hold off;


figure;
errorbar(2:K,db(:,1),db(:,2)-db(:,1),db(:,1)-db(:,3));
xlabel('k');
ylabel('Davies-Buldin Index');
xlim([2,K]);
hold on;
plot([2,K],[1,1],'k--');
hold off;

[qx,qy] = max(Q(logical(db(:,1)<1)));
figure;hold on;
plot(2:K,(Q-mean(Q))/(max(Q)-min(Q)),'-o','MarkerSize',6,'MarkerFaceColor','b');
plot(2:K,(db(:,1)-mean(db(:,1)))/(max(db(:,1))-min(db(:,1))),'-^','MarkerSize',6,'MarkerFaceColor','r');
xlabel('k');
ylabel('Modularity');

plot(qy+1,qx,'ks','MarkerSize',12)
hold off;

d = diag(DA);
[dx,dy]=sort(abs(d),'descend');
dx
dy
K
size(d)
figure;plot(1:K,d(dy),'-o','MarkerSize',6,'MarkerFaceColor','b');
xlabel('k');
ylabel('Eigenvalues');
hold on;
plot([1,K],[0,0],'k--');
hold off;



end

