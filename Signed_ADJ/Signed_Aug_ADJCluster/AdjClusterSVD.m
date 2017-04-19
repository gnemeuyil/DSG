function [Jf1, kf1, db1, Q1, ang1, acc1,Cf] = AdjClusterSVD(A,K,L,c,Partition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AdjCluster(A) or AdjCluster(A,K,L,c) or AdjCluster(A,K,L,c,Partition)% 
%Default K = min(100, ceil(m/10)); L = 2; c = 0; p = 'on';            %
%Input:                                                               %
%   A: Adjacency matrix                                               %
%   K: Maximal number of eigenpair to try                             %
%       Default: minimal between 100 and 1/50 of A's rank             %
%   L: K-mean algorithm's parameter distance                          %
%   c: min acceptable distance to the original                        %
%   p: plot option                                                    %
%   Partition: Known partition information to calculate the accuracy  %
%                                                                     %
%Output:                                                              %
%   Jf: Partition by AdjCluster                                       %
%   Cf: Centroids of each cluster in the spectral space               %
%   kf: best choice of k by different measures                        %
%       kf = [kf1 kf2 kf3]                                            %
%       kf1: best choice of k by difference of edges in and outside   %  
%       the communities                                               %
%       kf2: best choise of k by DB index                             %
%       kf3: best choise of k by modularity                           %
%   db: DB index                                                      %
%   Q: Modularity                                                     %
%   acc: Compare the known partition with partition by AdjCluster     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = size(A,1);
if ~exist('K','var') && ~exist('L','var')&& ~exist('c','var')
    K = min(100, ceil(m/10));
    L = 2;
    c = 0;
end;
    p = 'on';
    p='off';

% [V,D] = eigs(A,2,'lm');
% if D(1,1)*D(2,2)<0
%     s = 'lm';
% else s = 'la';
% end;
tic
%[V,D] = eigs(A,K,'la');
%d = nnz(sum(D)>0);
d=K;

db = zeros(min(d,K)-1,3);
kf = zeros(3,1);
Q = zeros(min(d,K)-1,1);
ang = zeros(min(d,K)-1,3);
acc = zeros(min(d,K)-1,1);
temp = [10 0];

for i = 2:min(d,K)
    %[V,D]=eigs(A,i,'lm');
    [U,S,V]=svds(A,i,'L');
    [Jt, dbt, ang3, Ct] = AdjCl(V(:,1:i),L,c);
    
    db(i-1,:) = dbt;
    ang(i-1,:) = ang3;
    if db(i-1,1)<temp(1)
        temp(1) = db(i-1,1);
        Jf{1} = Jt;
        kf(1) = i;
        Cf{1} = Ct;
    elseif temp(1) == db(i-1,1)
        WARNING('Multiple Solution!)');
    end;
    
    
    
    P = idx2lgc(Jt);
    Q(i-1) = SignQfunction(A,P);
%     if strcmp(s, 'lm') && min(min(A))==0
%         Q(i-1)=-Q(i-1);
%     end;
    if Q(i-1)>temp(2)
        temp(2) = Q(i-1);
        Jf{2} = Jt;
        kf(2) = i;
        Cf{2} = Ct;
    elseif temp(2) == Q(i-1)
        WARNING('Multiple Solution!)');
    end;
    
    if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
        acc(i-1) = PartitionAccuracy(P, Partition);
    end
end;
toc;
t=toc;


if strcmp(p,'on')
    figure;plot(2:min(d,K),db(:,1)','-o','MarkerSize',6,'MarkerFaceColor','b');
    xlabel('k','FontSize',14);
    ylabel('Davies-Buldin Index','FontSize',14);
    %hold on;
    %plot([2,K],[1,1],'k--');
    %hold off;
    box on;
    axis square;
    
%     figure;errorbar(2:K,db(:,1),db(:,2)-db(:,1),db(:,1)-db(:,3));
%     xlabel('k');
%     ylabel('Davies-Buldin Index');
%     xlim([2,K]);
%     hold on;
%     plot([2,K],[1,1],'k--');
%     hold off;
%     box on;
%     axis square;
% 
%     figure;errorbar(2:K,ang(:,2),ang(:,2)-db(:,1),db(:,3)-db(:,2));
%     xlabel('k');
%     ylabel('Angles between to the centroids');
%     xlim([2,K]);
%     hold on;
%     plot([2,K],[1,1],'k--');
%     hold off;
%     box on;
%     axis square;
% 
%     [qx,qy] = max(Q(logical(db(:,1)<1)));
    figure;plot(2:K,abs(Q),'-o','MarkerSize',6,'MarkerFaceColor','b');
    xlabel('k','FontSize',14);
    ylabel('Modularity','FontSize',14);
%     hold on;
%     plot(qy+1,qx,'ks','MarkerSize',12)
%     hold off;
    box on;
    axis square;

    % d = diag(D);
    % [dx,dy]=sort(abs(d),'descend');
%     figure;plot(1:K,diag(D),'-o','MarkerSize',6,'MarkerFaceColor','b');
%     xlabel('k');
%     ylabel('Eigenvalues');
%     if min(diag(D))<0
%         hold on;
%         plot([1,K],[0,0],'k--');
%         hold off;
%     end;
%     xlim([1,K]);
%     box on;
%     axis square;
end;

Jf1 = Jf;
kf1 = kf;
db1 = db;
Q1 = Q;
ang1 = ang;
acc1 = acc;


%[V,D] = eigs(A,K,'lm');
[U,S,V]=svds(A,K,'L');

d = nnz(sum(S)>0);
if d >K%close this loop
db = zeros(min(K-d,K-1),3);
kf = zeros(3,1);
Q = zeros(min(K-d,K-1),1);
ang = zeros(min(K-d,K-1),3);
acc = zeros(min(K-d,K-1),1);
temp = [10 0];

for i = 1:min(K-d,K-1)
    %[V,D]=eigs(A,i,'lm');
    [Jt, dbt, ang3] = AdjCl(V(:,[1,(K-i+1):(K)]),L,c);
    %[Jt, dbt, ang3] = AdjCl1(V(:,(K-i+1):(K)),L,c);
    
    db(i,:) = dbt;
    ang(i,:) = ang3;
    if db(i,1)<temp(1)
        temp(1) = db(i,1);
        Jf{1} = Jt;
        kf(1) = i;
    elseif temp(1) == db(i,1)
        WARNING('Multiple Solution!');
    end;
    
    P = idx2lgc(Jt);
    Q(i) = SignQfunction(A,P);
%     if strcmp(s, 'lm') && min(min(A))==0
%         Q(i-1)=-Q(i-1);
%     end;
    if Q(i)>temp(2)
        temp(2) = Q(i);
        Jf{2} = Jt;
        kf(2) = i;
    elseif temp(2) == Q(i)
        WARNING('Multiple Solution!)');
    end;
    
    if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
        acc(i) = PartitionAccuracy(P, Partition);
        size(P)
    end
end;


if strcmp(p,'on')
    figure;plot(2:(min(K-d,K-1)+1),db(:,1)','-o','MarkerSize',6,'MarkerFaceColor','b');
    xlabel('k');
    ylabel('Davies-Buldin Index');
    %hold on;
    %plot([2,K],[1,1],'k--');
    %hold off;
    box on;
    axis square;
    
%     figure;errorbar(2:K,db(:,1),db(:,2)-db(:,1),db(:,1)-db(:,3));
%     xlabel('k');
%     ylabel('Davies-Buldin Index');
%     xlim([2,K]);
%     hold on;
%     plot([2,K],[1,1],'k--');
%     hold off;
%     box on;
%     axis square;
% 
%     figure;errorbar(2:K,ang(:,2),ang(:,2)-db(:,1),db(:,3)-db(:,2));
%     xlabel('k');
%     ylabel('Angles between to the centroids');
%     xlim([2,K]);
%     hold on;
%     plot([2,K],[1,1],'k--');
%     hold off;
%     box on;
%     axis square;
% 
%     [qx,qy] = max(Q(logical(db(:,1)<1)));
%     figure;plot(2:K,Q,'-o','MarkerSize',6,'MarkerFaceColor','b');
%     xlabel('k');
%     ylabel('Modularity');
%     hold on;
%     plot(qy+1,qx,'ks','MarkerSize',12)
%     hold off;
%     box on;
%     axis square;

    % d = diag(D);
    % [dx,dy]=sort(abs(d),'descend');
%     figure;plot(1:K,diag(D),'-o','MarkerSize',6,'MarkerFaceColor','b');
%     xlabel('k');
%     ylabel('Eigenvalues');
%     if min(diag(D))<0
%         hold on;
%         plot([1,K],[0,0],'k--');
%         hold off;
%     end;
%     xlim([1,K]);
%     box on;
%     axis square;
end;

Jf2 = Jf;
kf2 = kf;
db2 = db;
Q2 = Q;
ang2 = ang;
acc2 = acc;

end;





