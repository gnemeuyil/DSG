k0 = [200,180,170,150,140];
% k = [200000,180000,170000,150000,140000];
k=k*100;
k=k*250;
k=k*500;
% signed graph with xmin =5 and d(1)=alpha=2.3, d(2)=0.15
 
tic
k=k0*300;
[A B]=Bal_N(k,2,[2.3,0.15]);
toc
pos=nnz(A);
neg=nnz(B);
 
k1=k0;
tic
[A1 B1]=Bal_N(k1,2,[2.3,0.12]);
toc
pos=nnz(A1)
neg=nnz(B1)
 
 
k1=k0;
tic
[A2 B2]=Bal_N(k1,5,[2.3,2.3]);
toc
pos=nnz(A2>0)
neg=nnz(B2<0)
%generate random graph with unclear community structure
[A B]=Bal_N(k,1,[0.3,0.25]);
M=-B;
A1 = A+M;
M = triu(M,1);
nn=find(M==1);
y=randsample(nn,ceil(0.50*size(nn,1)));
M0 = M;
M0(y)=-M0(y);
M0 = M0+M0';
A2 = A+M0;
%generate a random k partite graph
p=[0,0.1,0.1,0.1,0.1;
    0,0,0.1,0.1,0.1;
    0,0,0,0.1,0.1;
    0,0,0,0,0.1;
    0,0,0,0,0]*1+triu(rand(5)*0.2,1)
tic;[Ap, partition, q] = KpartiteGenerator(k,p);toc;q
tic;[Jf0 kf0 db0 Q0 ang0 acc0] = AdjCluster(Ap,10,2,0,partition);toc;
%go to the UniAdjcluster_Package folder
[Jf1 kf1 db1 Q1 ang1 acc1] = AdjCluster(A1,10,2,0,partition);
[Jf2 kf2 db2 Q2 ang2 acc2] = AdjCluster(A2,10,2,0,partition);
tic;[Jf0 kf0 db0 Q0 ang0 acc0] = AdjCluster(Ap,10,2,0,partition);toc;
 
 
%%
alpha=2.3;
xmin=5;
n=100;
A=rand(n);
d = round(xmin*(1-rand(n,1)).^(1/(alpha-1)));
P = d*d'/sum(d);
R = logical(A<=P);
R = triu(R,1);
R = R + R';
nnz(R)
 
 
%% PLUS2
% 0.2 negative edges
k0 = [200,180,170,150,140];
alpha=1.3;
d=[alpha,alpha];
P=0.2;
k=k0*600;
[A B]=Bal_N(k,5,d,P);
% partition
kn = max(size(k));%# of comm
kcum = [0,cumsum(k)];
n = kcum(kn+1);
partition = sparse(n,kn);
for i =1:kn
    partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
end;
partition = logical(partition);
tic;[Jf25 kf25 db25 Q25 ang25 acc25] = AdjCluster(A+B,10,2,0,partition);toc;
tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+B,10,2,0,partition);toc
% 0.5 negative edges
k0 = [200,180,170,150,140];
%alpha=1.3;
d=[1.3,1.3];
P=0.2;
k=k0*1200;%*1195;
[A, B]=Bal_N(k,5,d,P);
% partition
kn = max(size(k));%# of comm
kcum = [0,cumsum(k)];
n = kcum(kn+1);
partition = sparse(n,kn);
for i =1:kn
    partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
end;
partition = logical(partition);
tic;[Jf20, kf20, db20, Q20, ang20, acc20] = AdjCluster(A+B,10,2,0,partition);toc;
tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+B,10,2,0,partition);toc;
acc
acc20
 
[ Acc1 ] = Accuracy( Jf1,Jf,2);
 
nnz(B)/(sum(k)^2-k*k')
nnz(A)/(k*k')
sum(sum(B>0))
sum(sum(B<0))
 
Close all;
%% K-partite
clear all
k0 = [200,180,170,150,140];
k=k0*600;
%p=(triu(ones(5,5)*0.1,1)+triu(rand(5)*0.2,1)).*triu(ones(5,1)*(10./k),1);
p=(triu(ones(5,5)*0.2,1)).*triu(ones(5,1)*(10./k),1);
p=p+p';
tic;[Ap,A partition, q] = KpartiteGenerator(k,p);toc;
%tic;[Ap,A partition] = KpartiteGenerator2(k);toc;
tic;[Jf0 kf0 db0 Q0 ang0 acc0] = AdjCluster(Ap+A,10,2,0,partition);toc;
tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+Ap,10,2,0,partition);toc;
 
Q0
acc0
db0
ang0
Q
acc
 
%% Epinion
X=X+1;
Y=Y+1;
IN=E<0;
IP=E>0;
clear E
n=max(max(X),max(Y));
A=sparse(X(IP),Y(IP),1,n,n);
A=double(A|A');
N=sparse(X(IN),Y(IN),-1,n,n);
clear X Y IN IP 
tic;A(N|N')=-1;toc; %201
[x,y]=find(diag(A));
A=A-sparse(1:n,1:n,diag(A),n,n);
clear N
tic;[Jf0 kf0 db0 Q0 ang0 acc0] = AdjCluster(A,50,2,0);toc;
 
% find detailed information about clusters 
I1=(Jf0{2}==1);
I2=(Jf0{2}==2);
I3=(Jf0{2}==3);
I4=(Jf0{2}==4);
I5=(Jf0{2}==5);
nnz(I1)
nnz(I2)
nnz(I3)
nnz(I4)
nnz(I5)
 
D1=A(I1,I1);
nnz(D1)
D2=A(I2,I2);
nnz(D2)
D3=A(I3,I3);
nnz(D3)
D4=A(I4,I4);
nnz(D4)
D5=A(I5,I5);
nnz(D5)
 
 
%% input files for Python
[X,Y,E]=find(B+A);
[Xpar,Ypar]=find(partition);
dlmwrite('PL_edge_1000k_20.txt', [X Y E],'delimiter', ' ','precision','%10.10g');
dlmwrite('PL_partition_1000k_20.txt', [Xpar Ypar],'delimiter', ' ','precision','%10.10g');
 
 
%% compute accuracy
L1 = idx2lgc(VarName1);
[n,m]=size(L1);
L2=sparse(n,m);
L2(:,1:max(VarName2))=idx2lgc(VarName2);%bench mark
acc = PartitionAccuracy(L1, L2);
 
L2=idx2lgc(VarName2);%bench mark
[n,m]=size(L2);
L0 = idx2lgc(VarName1);
L1=sparse(n,m);
L1(:,1:4)=L0(:,1:4);
L1(:,5)=sum(L0(:,5:size(L1,2)),2);
acc = PartitionAccuracy(L1, L2);
 
 
 
 
 
%% Signed AUG_ADJ for Asymmetric Adjacency Matrices
clear
% 0.2 negative edges
k0 = [240,220,200,180,160];
%alpha=1.3;
d=[1.3,1.3];
P=0.5;
k=k0*1000;
 
T=zeros(50,3);
ACC=zeros(50,3);
S=zeros(50,4);
ANG=zeros(50,3);
DBI=zeros(50,3);
MOD=zeros(50,1);
for j=1:50
 
    [A, B]=Bal_N(k,5,d,P);
    % partition
    kn = max(size(k));%# of comm
    kcum = [0,cumsum(k)];
    n = kcum(kn+1);
    partition = sparse(n,kn);
    for i =1:kn
        partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
    end;
    partition = logical(partition);
    tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Generalized_ADJCluster(A+B,10,0,partition,1);toc;
    t1=toc;
    tic;[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(A+B,10,2,0,partition);toc;
    t2=toc;
    tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+B,10,2,0,partition);toc;
    t3=toc;
    T(j,:)=[t1 t2 t3];
    ACC(j,:)=[max(acc1) acc0(4) acc(4)];
    S(j,:)=[nnz((A+B)>0) nnz((A+B)<0) nnz(A) nnz(B)];
    ANG(j,:)=ang1(5,:);
    DBI(j,:)=db1(5,:);
    MOD(j,:)=Q1(5);
end;
 
dlmwrite('S.txt', S,'delimiter', ' ','precision','%10.10g');
dlmwrite('ACC.txt', ACC,'delimiter', ' ','precision','%10.10g');
dlmwrite('ANG.txt', ANG,'delimiter', ' ','precision','%10.10g');
dlmwrite('DBI.txt', DBI,'delimiter', ' ','precision','%10.10g');
dlmwrite('MOD.txt', MOD,'delimiter', ' ','precision','%10.10g');
dlmwrite('T.txt', T,'delimiter', ' ','precision','%10.10g');
%% Accuracy and run time comparison
R=T;
 
figure;
bar(1:size(R,1), R, 1)
% Set the axis limits
axis([0 size(R,1)+1 0 max(max(R))+2])
set(gca, 'XTick', 1:size(R,1)+1)
 
 
figure;bar(1:3,mean(R),0.5)
hold on;
h=errorbar(1:3,mean(R),std(R),'c'); set(h,'linestyle','none')
 
%% (Directed)Reconstructed Demostration Graph with negative inter community edge
X=[1 1 1 2 2 3 3 3 4 4 5 5 6 6 7 8 8 9 9 9 10 11 11 11 12 13 14 14 15 16 16 16 16 17 17 17 18 19 20 20 21 21 22 23 23 23 24 25 25 25 25]
Y=[2,5,6,4,5,4,7,8,1,14,8,15,3,7,5,14,15,10,12,13,12,9,12,13,13,10,1,7,2,18,20,22,23,19,20,24,21,18,19,22,24,25,17,18,21,24,20,16,18,19,20]
JF=[1;1;1;1;1;1;1;1;3;3;3;3;3;1;1;2;2;2;2;2;2;2;2;2;2]
M{1}='g*';
M{2}='m*';
 
G=sparse(X,Y,1,max(X),max(X));
[v,d]=eigs(G,3,'LR');
 
% compare positive and negative intercluster
% edge 5 to 10 is added from cluster 1 to 3, so the eigenpair 
% corresponding to cluster 3 is perturbed, K=3
% v3_5_10 means the 3rd eigen vector of v_5_10
 
% 5+>10
[v3_5_10,d3_5_10,v_5_10,d_5_10, G_5_10 ] = getEigPair(G,5,10,1,3,3);
plotClusters( JF,v_5_10 ,[1 3],[5,10],M )
 
% 5->10
[v3_5_10_neg,d3_5_10_neg,v_5_10_neg,d_5_10_neg, G_5_10_neg ] = getEigPair(G,5,10,-1,3,3);
plotClusters( JF,v_5_10_neg ,[1 3],[5,10],M )
 
% edge 10 to 5, so the eigenpair corresponding to cluster 1 is perturbed, K=1
% 10+>5
[v1_10_5,d1_10_5,v_10_5,d_10_5, G_10_5 ] = getEigPair(G,10,5,1,3,1);
plotClusters( JF,v_10_5 ,[1 3],[5,10],M )
 
% 10->5
[v1_10_5_neg,d1_10_5_neg,v_10_5_neg,d_10_5_neg, G_10_5_neg ] = getEigPair(G,10,5,-1,3,1);
plotClusters( JF,v_10_5_neg,[1 3],[5,10],M )
 
% adding edges 8 and 25
% 8+>25
[v3_8_25,d3_8_25,v_8_25,d_8_25, G_8_25 ] = getEigPair(G,8,25,1,3,2);
plotClusters( JF,v_8_25 ,[1 2],[8,25],M )
 
% 8->25
[v3_8_25_neg,d3_8_25_neg,v_8_25_neg,d_8_25_neg, G_8_25_neg ] = getEigPair(G,8,25,-1,3,2);
plotClusters( JF,v_8_25_neg ,[1 2],[8,25],M )
 
% 25+>8
[v3_25_8,d3_25_8,v_25_8,d_25_8, G_25_8 ] = getEigPair(G,25,8,1,3,2);
plotClusters( JF,v_25_8 ,[1 2],[8,25],M )
 
% 25->8
[v3_25_8_neg,d3_25_8_neg,v_25_8_neg,d_25_8_neg, G_25_8_neg ] = getEigPair(G,25,8,-1,3,2);
plotClusters( JF,v_25_8_neg ,[1 2],[8,25],M )
 
[v(:,2) v_8_25(:,2) v_8_25_neg(:,2)]
 
G258=G_25_8;
G258(8,25)=1;
[v258,d258]=eigs(G258,3,'lr')
plotClusters( JF,v258,[1 2],[8,25],M )
 
G258neg=G258;
G258neg(25,8)=-1;
G258neg(8,25)=-1;
[v258neg,d258neg]=eigs(G258neg,3,'lr')
plotClusters( JF,v258neg ,[1 2],[8,25],M )
 
[v3_8_25_neg,d3_8_25_neg,v_8_25_neg,d_8_25_neg, G_8_25_neg ] = getEigPair(G,8,25,-1,3,2);
plotClusters( JF,v_8_25_neg ,[1 2],[8,25],M )
[v3_25_8_neg,d3_25_8_neg,v_25_8_neg,d_25_8_neg, G_25_8_neg ] = getEigPair(G,25,8,-1,3,2);
plotClusters( JF,v_25_8_neg ,[1 2],[8,25],M )
[v258neg,d258neg]=eigs(G258neg,3,'lr');
plotClusters( JF,v258neg ,[1 2],[8,25],M )
 
 
% adding edges 8 and 20
% 8+>20
[v3_8_20,d3_8_20,v_8_20,d_8_20, G_8_20 ] = getEigPair(G,8,20,1,3,2);
plotClusters( JF,v_8_20 ,[1 2],[8,20],M )
 
% 8->20
[v3_8_20_neg,d3_8_20_neg,v_8_20_neg,d_8_20_neg, G_8_20_neg ] = getEigPair(G,8,20,-1,3,2);
plotClusters( JF,v_8_20_neg ,[1 2],[8,20],M )
 
 
 
G_25_8=G;
G_25_8(25,8)=1;
[v258,d258]=eigs(G_25_8,3,'LR')
plotClusters( JF,v258 ,[1 2],[],[] )
G_25_8_neg=G;
G_25_8_neg(25,8)=-1;
[v258neg,d258neg]=eigs(G_25_8_neg,3,'LR')
plotClusters( JF,v258neg ,[1 2],[],[] )
 
%comparison on 10->5 +/-
[v,d]=eigs(G,3,'LR')
[v105,d105]=eigs(G_10_5,3,'LR')
[v105neg,d105neg]=eigs(G_10_5_neg,3,'LR')
[v(:,1) v(:,3) v105(:,1) v105(:,3) v105neg(:,1) v105neg(:,3)]
 
%comparison on 25->8 +/-
[v,d]=eigs(G,3,'LR')
[v258,d258]=eigs(G_25_8,3,'LR')
[v258neg,d258neg]=eigs(G_25_8_neg,3,'LR')
[v(:,1) v(:,2) v258(:,1) v258(:,2) v258neg(:,1) v258neg(:,2)]
 
%plots for multiple intercluster negative edges
M{1}='g*';
M{2}='m*';
%M{3}='y*';
%M{4}='c*';
N=[10,13];
G_10_5_13_8_neg=G;
G_10_5_13_8_neg(10,5)=-1;
G_10_5_13_8_neg(13,8)=-1;
N=[10,13];
[v105138neg,d105138neg]=eigs(G_10_5_13_8_neg,3,'LR')
plotClusters( JF,v105138neg ,[1 3],N,M);
plotClusters( JF,v105138neg ,[1 2 3],N,M);
 
% negative inner cluster perturbations
G_group1_7_15_neg=G;
G_group1_7_15_neg(7,15)=-1;
[v_group1_7_15_neg,d_group1_7_15_neg]=eigs(G_group1_7_15_neg,3,'LR');
%plotClusters( JF,v_group1715_neg,[1 2 3],[],[] )
[v(:,1) v_group1_7_15_neg(:,1)]
 
G_group1_7_15=G;
G_group1_7_15(7,15)=1;
[v_group1_7_15,d_group1_7_15]=eigs(G_group1_7_15,3,'LR');
%plotClusters( JF,v_group1715_neg,[1 2 3],[],[] )
[v(:,1) v_group1_7_15(:,1) v_group1_7_15_neg(:,1)]
 
G_group1_7_15_neg2=G;
G_group1_7_15_neg2(7,15)=-2;
[v_group1_7_15_neg2,d_group1_7_15_neg2]=eigs(G_group1_7_15_neg2,3,'LR');
%plotClusters( JF,v_group1715_neg,[1 2 3],[],[] )
[v(:,1) v_group1_7_15_neg(:,1) v_group1_7_15_neg2(:,2)]
[d(1,1) d_group1_7_15_neg(1,1) d_group1_7_15_neg2(2,2)]
 
%% (Directed)Reconstructed Demostration Graph with negative INNER community edges to test COMPLEX eigenpairs and their properties
% add random negative edges to graph and check for complex eigenpair plots
% only for group 1
 
X=[1 1 1 2 2 3 3 3 4 4 5 5 6 6 7 8 8 9 9 9 10 11 11 11 12 13 14 14 15 16 16 16 16 17 17 17 18 19 20 20 21 21 22 23 23 23 24 25 25 25 25];
Y=[2,5,6,4,5,4,7,8,1,14,8,15,3,7,5,14,15,10,12,13,12,9,12,13,13,10,1,7,2,18,20,22,23,19,20,24,21,18,19,22,24,25,17,18,21,24,20,16,18,19,20];
JF=[1;1;1;1;1;1;1;1;3;3;3;3;3;1;1;2;2;2;2;2;2;2;2;2;2];
G=sparse(X,Y,1,max(X),max(X));
idx=[1:8, 14,15];
sub_G=G(idx,idx);
M{1}='g*';
M{2}='m*';
I=1:2;
[d,v]=eigs(G,10,'LR');
plotClusters( JF,[abs(d(:,1)) d(:,2)] ,I,N,M)

M{3}='y*';
M{4}='c*';
M{5}='b*';
 
% the sub graph of graoup 1 has only 20 edges, contains node 1:8, 14 and 15
% the size of the matrix is 10 by 10
subG=sub_G;
K=10;
val=-1;
[ D_subG_LR_tilde,V_subG_LR_tilde,D_subG_SR_tilde,V_subG_SR_tilde, subG_tilde ,S] = EigenPerturb( subG, K,val);
 
for i=1:K
    [D_subG_LR_tilde{i},D_subG_SR_tilde{i}]
end;
 
figure;
plot( V_subG_tilde{4}(:,1), V_subG_tilde{4}(:,6),'g*' );
 
 
%% Generating Dynamic graphs to analyze perturbations of negative edges
X=[1 1 1 2 2 3 3 3 4 4 5 5 6 6 7 8 8 9 9 10 11 11 11 12 13 13 13 14 15 16 16 16 16 17 17 17 18 19 20 20 21 21 22 23 23 23 24 25 25 25 25]
Y=[2,5,6,4,5,4,7,8,1,9,8,10,3,7,5,9,10,1,7,2,12,14,15,14,11,14,15,15,12,18,20,22,23,19,20,24,21,18,19,22,24,25,17,18,21,24,20,16,18,19,20]
JF=[1;1;1;1;1;1;1;1;1;1;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2]
G=sparse(X,Y,1,max(X),max(X));
idx1=1:10;
idx2=11:15;
idx3=16:25;
sub_G1=G(idx1,idx1);
sub_G2=G(idx2,idx2);
sub_G3=G(idx3,idx3);
M{1}='g*';
M{2}='m*';
M{3}='y*';
M{4}='c*';
M{5}='b*';
subG=sub_G1;
K=10;
val=-1;%could be negative or positive
[ D_subG_LR_tilde,V_subG_LR_tilde,D_subG_SR_tilde,V_subG_SR_tilde, subG_tilde ,S] = EigenPerturb( subG, K,val);
StartDate=100;
[ T_Nodes, T_Edges ] = SavetoFileGephi( subG_tilde,StartDate,S,'Node.txt', 'Edge.txt');
for i=1:K
    [D_subG_LR_tilde{i},abs(D_subG_LR_tilde{i}),D_subG_SR_tilde{i},abs(D_subG_SR_tilde{i})]
end;
 
%testing for existance of Perrron root
for i=1:50
   [ D_subG_LR_tilde,V_subG_LR_tilde,D_subG_SR_tilde,V_subG_SR_tilde, subG_tilde ,S] = EigenPerturb( subG, K,val);
   RealEigVal(i)=D_subG_LR_tilde{K}(1)==real(D_subG_LR_tilde{K}(1));
   %RealEigVec(i)=sum(V_subG_LR_tilde{K}(:,1)==real(V_subG_LR_tilde{K}(:,1)));
   PosEigVec(i)=(sum(V_subG_LR_tilde{K}(:,1)==abs(V_subG_LR_tilde{K}(:,1)))==10)|(sum(V_subG_LR_tilde{K}(:,1)==-abs(V_subG_LR_tilde{K}(:,1)))==10);
   SS(i)=S;
   [D_subG_LR_tilde,V_subG_LR_tilde,D_subG_SR_tilde,V_subG_SR_tilde]
end;
[RealEigVal' PosEigVec']
 
%% Test Moduli and Split methods
X=[1 1 1 2 2 3 3 3 4 4 5 5 6 6 7 8 8 9 9 10 11 11 11 12 13 13 13 14 15 16 16 16 16 17 17 17 18 19 20 20 21 21 22 23 23 23 24 25 25 25 25]
Y=[2,5,6,4,5,4,7,8,1,9,8,10,3,7,5,9,10,1,7,2,12,14,15,14,11,14,15,15,12,18,20,22,23,19,20,24,21,18,19,22,24,25,17,18,21,24,20,16,18,19,20]
JF=[1;1;1;1;1;1;1;1;1;1;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2]
G=sparse(X,Y,1,max(X),max(X));
idx1=1:10;
idx2=11:15;
idx3=16:25;
sub_G1=G(idx1,idx1);
sub_G2=G(idx2,idx2);
sub_G3=G(idx3,idx3);
M{1}='g*';
M{2}='m*';
M{3}='y*';
M{4}='c*';
M{5}='b*';
K=1;
clear XY
val=-1;%could be negative or positive
Partition=idx2lgc(JF);
for i=1:100
[ ~,~,~,~, subG_tilde1 ,S1] = EigenPerturb( sub_G1, K,val);
[ ~,~,~,~,  subG_tilde2 ,S2] = EigenPerturb( sub_G2, K,val);
[ ~,~,~,~,  subG_tilde3 ,S3] = EigenPerturb( sub_G3, K,val);
 
G_tilde=sparse(25,25,0);
 
G_tilde(1:10,1:10)=subG_tilde1;
G_tilde(11:15,11:15)=subG_tilde2;
G_tilde(16:25,16:25)=subG_tilde3;
[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(G_tilde,23,0,1,0,Partition);%split
[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(G_tilde,23,0,1,1,Partition);%mudulus
 
ACC(i,:)=[max(acc1) max(acc2)];
% [Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Generalized_ADJCluster(G_tilde,23,0,1,0,3);%split
% [Jf4, kf4, db4, Q4, ang4, acc4,Cf4,v42,d42,v42,d42,I4,cand4] = Generalized_ADJCluster(G_tilde,23,0,1,1,3);%mudulus
JFS(:,[2*i-1,2*i])=[Jf1{2} Jf2{2}];
 
%[x,y]=find(G_tilde<0);
%XY(:,[2*i-1,2*i])=[x,y];
end
 
%% Test Complex Eigenvector vs Graph Structure
X=[1 1 1 2 2 3 3 3 4 4 5 5 6 6 7 8 8 9 9 10 11 11 11 12 13 13 13 14 15 16 16 16 16 17 17 17 18 19 20 20 21 21 22 23 23 23 24 25 25 25 25]
Y=[2,5,6,4,5,4,7,8,1,9,8,10,3,7,5,9,10,1,7,2,12,14,15,14,11,14,15,15,12,18,20,22,23,19,20,24,21,18,19,22,24,25,17,18,21,24,20,16,18,19,20]
JF=[1;1;1;1;1;1;1;1;1;1;3;3;3;3;3;2;2;2;2;2;2;2;2;2;2]
G=sparse(X,Y,1,max(X),max(X));
idx1=1:10;
idx2=11:15;
idx3=16:25;
sub_G1=G(idx1,idx1);
sub_G2=G(idx2,idx2);
sub_G3=G(idx3,idx3);
M{1}='g*';
M{2}='m*';
M{3}='y*';
M{4}='c*';
M{5}='b*';
 
clear X Y
 
val=-1;%could be negative or positive
Partition=idx2lgc(JF);
K=5;
T=500;
Theta=10^(-8);
SS=[];
R=zeros(T,K);
%R1=zeros(T,1);
%testing for existance of Perrron root
for i=1:T
    [ D_subG_LR_tilde,V_subG_LR_tilde,D_subG_SR_tilde,V_subG_SR_tilde, subG_tilde ,S] = EigenPerturb( subG, K,val);
    RealEigVal(i)=D_subG_LR_tilde{K}(1)==real(D_subG_LR_tilde{K}(1));
    %RealEigVec(i)=sum(V_subG_LR_tilde{K}(:,1)==real(V_subG_LR_tilde{K}(:,1)));
    PosEigVec(i)=(sum(V_subG_LR_tilde{K}(:,1)==abs(V_subG_LR_tilde{K}(:,1)))==10)|(sum(V_subG_LR_tilde{K}(:,1)==-abs(V_subG_LR_tilde{K}(:,1)))==10);
    SS=[SS S'];
    [I,J]=ind2sub(size(subG_tilde),S)
    for j=1:K    
        [~,N1]=max(abs(D_subG_LR_tilde{j}));
        [~,N]=max(abs(D_subG_SR_tilde{j}));
        [D_subG_LR_tilde{j},abs(D_subG_LR_tilde{j}),D_subG_SR_tilde{j},abs(D_subG_SR_tilde{j})]
        warning('Eigenvectors corresponding to spectral radii')
        [V_subG_LR_tilde{j}(J,N1),abs(V_subG_LR_tilde{j}(J,N1)),V_subG_SR_tilde{j}(J,N),abs(V_subG_SR_tilde{j}(J,N))]% nodes with -> negative edge
        [V_subG_LR_tilde{j}(:,N1),abs(V_subG_LR_tilde{j}(:,N1)),V_subG_SR_tilde{j}(:,N),abs(V_subG_SR_tilde{j}(:,N))]% all the nodes
        %if nnz(D_subG_LR_tilde{j}==real(D_subG_LR_tilde{j}))>0&&nnz(D_subG_LR_tilde{j}(D_subG_LR_tilde{j}==real(D_subG_LR_tilde{j}))~=1)<2
        if nnz(D_subG_LR_tilde{j}==real(D_subG_LR_tilde{j}))>0&&nnz(D_subG_LR_tilde{j}(D_subG_LR_tilde{j}==real(D_subG_LR_tilde{j}))~=1)>=1
            R(i,j)=1;
        end;
%         if nnz(abs(D_subG_SR_tilde{j}-real(D_subG_SR_tilde{j}))<Theta)>0
%             R1(i)=R1(i)+1;
%         end;    
    end;
end;
nnz(R==K)
%nnz(R1==K)
sum(R)
 
%Percentage of the existence of at least one real positive eigenpair after each perturbation
y=sum(R,1)/10;
x=1:size(y,2);
figure;
hold on
bar(y,'BarWidth',0.4,'FaceColor',[0.5 0.5 0.5]);
ylim([0 max(y)+5]);
set(gca,'XTick',1:1:size(y,2));
for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
hold off
 
%Percentage distribution of the number of times of the existence of real positive eigenpair over all perturbations
y=[nnz(sum(R,2)==1) nnz(sum(R,2)==2) nnz(sum(R,2)==3) nnz(sum(R,2)==4) nnz(sum(R,2)==5)]./10;
x=1:size(y,2);
hold on
bar(y,'BarWidth',0.4,'FaceColor',[0.5 0.5 0.5]);
ylim([0 max(y)+5]);
set(gca,'XTick',1:1:size(y,2));
for i1=1:numel(y)
    text(x(i1),y(i1),num2str(y(i1),'%0.1f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end
hold off
 
% get the indices of pertrubation which caused the real eigenpair to
% disapear
[I,J]=ind2sub(size(subG_tilde),SS(:,R(:,2)==0))
[I,J]=ind2sub(size(subG_tilde),SS(:,R(:,3)==0))
 
%% Epinion Directed Signed
X=X+1;
Y=Y+1;
IN=E<0;
IP=E>0;
n=max(max(X),max(Y));
% 1 original matrix
A=sparse(X,Y,E,n,n);
 
 
 
% 2 use co-citation matrix
B=A*A';
B=B-sparse(1:n,1:n,1,n,n);
[Jf, kf, db, Q, ang, acc,Cf,v1,d1,v2,d2,I,cand] = Generalized_ADJCluster(B,30,0,[],1);
 
% find detailed information about clusters 
Num=[];
Density=[];
NumEdge=[];
NumEdgePos=[];
NumEdgeNeg=[];
for i=1:max(Jf{2})
    I=(Jf{2}==i);
    Num=[Num nnz(I)];
    D=A(I,I);
    NumEdge=[NumEdge nnz(D)];
    Density=[Density nnz(D)/nnz(I)^2];
    NumEdgePos=[NumEdgePos nnz(D>0)];
    NumEdgeNeg=[NumEdgeNeg nnz(D<0)];
end;
 
[Num;NumEdge;NumEdgePos;NumEdgeNeg]
Density
 
plotClusters( Jf{2},v1 ,[1 2 3],[],[] );
plotClusters( Jf{2},v1 ,[3 4 5],[],[] );
% 3 Skew Symmetric
 Z=(A-A')/2;
 [U,S,V]=svds(Z);
 US=U*S^(0.5);
 I=[];
%the plot takes a very~~~ long time
 PlotSkew2D(US,sqrt(S(1,1)/size(Z,1)),I,1);
 
 
 
 


