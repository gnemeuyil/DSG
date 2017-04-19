%% Signed AUG_ADJ for Asymmetric Adjacency Matrices With PFn Cluster
clear
% 0.2 negative edges
k0 = [240,220,200,180,160];
k=k0;
%k0 = [220,210,200,190,180];
%Inner cluster density for UNI is 0.2, Inter cluter for PL has alpha=1.3;
d=[0.4,0.2];
P=[0, 0.5];
rpt=1;
 
T=zeros(rpt,4);
ACC=zeros(rpt,4);
S=zeros(rpt,4);
ANG=zeros(rpt,3);
DBI=zeros(rpt,3);
MOD=zeros(rpt,1);
ANG2=zeros(rpt,3);
DBI2=zeros(rpt,3);
MOD2=zeros(rpt,1);
for j=1:rpt
 
    [A, B]=Bal_N(k,8,d,P);
    % partition
    kn = max(size(k));%# of comm
    kcum = [0,cumsum(k)];
    n = kcum(kn+1);
    partition = sparse(n,kn);
    for i =1:kn
        partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
    end;
    partition = logical(partition);
    tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(A+B,10,0,1,0,partition);toc;%split
    t1=toc;
    tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(A+B,10,0,1,1,partition);toc;%mudulus
    t2=toc;
    tic;[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(A+B,10,2,0,partition);toc;
    t3=toc;
    tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+B,10,2,0,partition);toc;
    t4=toc;
    T(j,:)=[t1 t2 t3 t4];
    ACC(j,:)=[max(acc1) max(acc2) acc0(4) acc(4)];
    S(j,:)=[nnz((A-B)>0) nnz((A-B)<0) nnz(A) nnz(B)];
    ANG(j,:)=ang1(5,:);
    DBI(j,:)=db1(5,:);
    MOD(j,:)=Q1(5);
    ANG2(j,:)=ang2(5,:);
    DBI2(j,:)=db2(5,:);
    MOD2(j,:)=Q2(5);
end;
ACC
 
 
dlmwrite('S.txt', S,'delimiter', ' ','precision','%10.10g');
dlmwrite('ACC.txt', ACC,'delimiter', ' ','precision','%10.10g');
dlmwrite('ANG.txt', ANG,'delimiter', ' ','precision','%10.10g');
dlmwrite('DBI.txt', DBI,'delimiter', ' ','precision','%10.10g');
dlmwrite('MOD.txt', MOD,'delimiter', ' ','precision','%10.10g');
dlmwrite('T.txt', T,'delimiter', ' ','precision','%10.10g');
 
%% Signed AUG_ADJ for DSG With non-PFn Cluster
clear
% 0.2 negative edges
k0 = [240,220,200,180,160];
k=k0;
%k0 = [220,210,200,190,180];
%Inner cluster density for UNI is 0.2, Inter cluter for UNI ;
d=[0.4,0.2];
P=[0, 1];
rpt=1;
 
T=zeros(rpt,4);
ACC=zeros(rpt,4);
S=zeros(rpt,4);
ANG=zeros(rpt,3);
DBI=zeros(rpt,3);
MOD=zeros(rpt,1);
ANG2=zeros(rpt,3);
DBI2=zeros(rpt,3);
MOD2=zeros(rpt,1);
for j=1:rpt
 
    [A, B]=Bal_N(k,9,d,P);
    % partition
    kn = max(size(k));%# of comm
    kcum = [0,cumsum(k)];
    n = kcum(kn+1);
    partition = sparse(n,kn);
    for i =1:kn
        partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
    end;
    partition = logical(partition);
    tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(A+B,20,0,1,0,partition);toc;%split
    t1=toc;
    tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(A+B,20,0,1,1,partition);toc;%mudulus
    t2=toc;
    tic;[Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Generalized_ADJCluster(A+B,20,0,1,2,partition);toc;%real
    t5=toc;
    tic;[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(A+B,20,2,0,partition);toc;
    t3=toc;
    tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+B,20,2,0,partition);toc;
    t4=toc;
    T(j,:)=[t1 t2 t3 t4 t5];
    ACC(j,:)=[max(acc1) max(acc2) acc0(4) acc(4)];
    S(j,:)=[nnz((A-B)>0) nnz((A-B)<0) nnz(A) nnz(B)];
    ANG(j,:)=ang1(5,:);
    DBI(j,:)=db1(5,:);
    MOD(j,:)=Q1(5);
    ANG2(j,:)=ang2(5,:);
    DBI2(j,:)=db2(5,:);
    MOD2(j,:)=Q2(5);
end;
ACC
[DBI(:,2) ANG(:,2) MOD]
T
[nnz(A) nnz(B)]
[max(Jf1{2}) max(Jf2{2}) max(Jf0{2}) max(Jf{2})]
 
%% Signed Tests for Aug_ADJ and Generalized_ADJ
clear
% 0.2 negative edges
k0 = [240,220,200,180,160];
k=k0;
%k0 = [220,210,200,190,180];
%Inner cluster density for UNI is 0.2, Inter cluter for PL has alpha=1.3;
d=[0.4,0.2];
P=[0.85, 0.5];
rpt=10;
 
% T=zeros(rpt,4);
ACC=zeros(rpt,5);
Clusters=zeros(rpt,5);
% S=zeros(rpt,4);
% ANG=zeros(rpt,3);
% DBI=zeros(rpt,3);
% MOD=zeros(rpt,1);
% ANG2=zeros(rpt,3);
% DBI2=zeros(rpt,3);
% MOD2=zeros(rpt,1);
for j=1:rpt
 
    [A, B]=Bal_N(k,9,d,P);
    % partition
    kn = max(size(k));%# of comm
    kcum = [0,cumsum(k)];
    n = kcum(kn+1);
    partition = sparse(n,kn);
    for i =1:kn
        partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
    end;
    partition = logical(partition);
    tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(A+B,20,0,1,0,partition);toc;%split
    t1=toc;
    tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(A+B,20,0,1,1,partition);toc;%mudulus
    t2=toc;
    tic;[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(A+B,20,2,0,partition);toc;
    t3=toc;
    tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(A+B,20,2,0,partition);toc;
    t4=toc;
    [Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Augmented_ADJCluster(A+B,20,1,0,0,1,partition);
       
%     T(j,:)=[t1 t2 t3 t4];
    ACC(j,:)=[max(acc1) max(acc2) max(acc3) acc0(max(Jf0{1}-1)) acc(max(Jf{1})-1)];
    Clusters(j,:)=[max(Jf1{2}) max(Jf2{2}) max(Jf3{2}) max(Jf0{1}) max(Jf{1})];
%     S(j,:)=[nnz((A-B)>0) nnz((A-B)<0) nnz(A) nnz(B)];
%     ANG(j,:)=ang1(5,:);
%     DBI(j,:)=db1(5,:);
%     MOD(j,:)=Q1(5);
%     ANG2(j,:)=ang2(5,:);
%     DBI2(j,:)=db2(5,:);
%     MOD2(j,:)=Q2(5);
end;
ACC
Clusters
%% Tests with GADJRe
clear
% 0.2 negative edges
k0 = [240,220,200,180,160];
k=k0;
%k0 = [220,210,200,190,180];
%Inner cluster density for UNI is 0.2, Inter cluter for UNI ;
d=[0.4,0.2];
P=[0.2, 1];
rpt=5;
T5=zeros(rpt,1);
T=zeros(rpt,4);
ACC=zeros(rpt,3);
S=zeros(rpt,4);
ANG=zeros(rpt,3);
DBI=zeros(rpt,3);
MOD=zeros(rpt,1);
ANG2=zeros(rpt,3);
DBI2=zeros(rpt,3);
MOD2=zeros(rpt,1);
 
for j=1:rpt
[A, B]=Bal_N(k,9,d,P);
    % partition
    kn = max(size(k));%# of comm
    kcum = [0,cumsum(k)];
    n = kcum(kn+1);
    partition = sparse(n,kn);
    for i =1:kn
        partition=partition+sparse((kcum(i)+1):(kcum(i+1)),i,1,n,kn);
    end;
    partition = logical(partition);
    tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(A+B,20,0,1,0,partition);toc;%split
    t1=toc;
    tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(A+B,20,0,1,1,partition);toc;%mudulus
    t2=toc;
    tic;[Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Generalized_ADJCluster(A+B,20,0,1,2,partition);toc;%real
    t5=toc;
    T5(j)=t5;
    ACC(j,:)=[max(acc1) max(acc2) max(acc3)];
end; 
 
 
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
 
%% Sampson
n=max([X;Y]);
% sort the ID and names
[~,SID]=sort(ID);
SortedNames=Names(SID);
 
Samp=sparse(X,Y,E.*(X~=Y),n,n);
 
[Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Augmented_ADJCluster(Samp,15,1,0,0,1);%AugAdj

tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Samp,15,0,1,0);toc;%split
t1=toc;
max(Jf1{2})
for i=1:max(Jf1{2})
    I=(Jf1{2}==i);
    SortedNames(I)   
end;
 
% comparsions
[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Samp,15,0,1,0);%split
max(Jf1{2})
Q1
 
IN=E<0;
IP=E>0;
n=max(max(X),max(Y));
SSamp=sparse(X(IP),Y(IP),1,n,n);
SSamp=double(SSamp|SSamp');
N=sparse(X(IN),Y(IN),-1,n,n);
clear IN IP 
tic;SSamp(N|N')=-1;toc; %201
[x,y]=find(diag(SSamp));
SSamp=SSamp-sparse(1:n,1:n,diag(SSamp),n,n);
clear N
 
[Jf, Cf, kf, db, Q, acc] = NormalCutS(SSamp,15,2,0);
 
[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(SSamp,15,2,0);
 
[max(Jf1{2}) max(Jf{1}) max(Jf0{1})]
max(Q1)
Q(max(Jf{1})-1)
Q0(max(Jf0{1})-1)
 
tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(Samp,15,0,1,1);toc;%mudulus
t2=toc
max(Jf2{2})
for i=1:max(Jf2{2})
    I=(Jf2{2}==i);
    SortedNames(I)   
end;
 
tic;[Jf, Cf, kf, db, Q, acc] = NormalCutS(Samp,15,2,0);toc;
max(Jf{1})
for i=1:max(Jf{1})
    I=(Jf{1}==i);
    SortedNames(I)   
end;
 
plotClustersComplex( Jf1{2},d11,v11 ,[1,2],[],[],SortedNames )
plotClustersComplex( Jf1{2},d11,v11 ,[1,3],[],[],SortedNames )
 
 
%% Epinion
X=X+1;
Y=Y+1;
n=max([X;Y]);
Epinion=sparse(X,Y,E.*(X~=Y),n,n);
[Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Augmented_ADJCluster(Epinion,50,1,0,0,1);%AugAdj
tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Epinion,50,0,1,0);toc;%split
t1=toc;
 
 
IN=0;
neg=0;
pos=0;
% detailed information of clustering results
for i=1:max(Jf1{2})
    I=(Jf1{2}==i);
    nnz(I)
    D=Epinion(I,I);
    nnz(D>0)
    nnz(D<0)   
    nnz(D)/(nnz(I)^2)
    IN=IN+nnz(D);
    pos=pos+nnz(D>0);
    neg=neg+nnz(D<0);
end;
IN
pos
neg
nnz(Epinion)-IN
nnz(Epinion>0)-pos
nnz(Epinion<0)-neg
max(Jf1{2})
 
 
 
[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Epinion,50,0,1,0);%split
max(Jf1{2})
Q1
 
IN=E<0;
IP=E>0;
n=max(max(X),max(Y));
SEpinion=sparse(X(IP),Y(IP),1,n,n);
SEpinion=double(SEpinion|SEpinion');
N=sparse(X(IN),Y(IN),-1,n,n);
clear IN IP 
tic;SEpinion(N|N')=-1;toc; %201
[x,y]=find(diag(SEpinion));
SEpinion=SEpinion-sparse(1:n,1:n,diag(SEpinion),n,n);
clear N
 
[Jf, Cf, kf, db, Q, acc] = NormalCutS(SEpinion,50,2,0);
 
[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(SEpinion,50,2,0);
 
[max(Jf1{2}) max(Jf{1}) max(Jf0{1})]
max(Q1)
Q(max(Jf{1})-1)
Q0(max(Jf0{1})-1)
 
 
 
 
 
 
%% Slashdot
n=max([X;Y]);
Slash=sparse(X,Y,E.*(X~=Y),n,n);
[Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Augmented_ADJCluster(Slash,50,1,0,0,1);%AugAdj
tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Slash,50,0,1,0);toc;%split
t1=toc;
max(Jf1{2})
 
IN=0;
neg=0;
pos=0;
% detailed information of clustering results
for i=1:max(Jf1{2})
    I=(Jf1{2}==i);
    nnz(I)
    D=Slash(I,I);
    nnz(D>0)
    nnz(D<0)   
    nnz(D)/(nnz(I)^2)
    IN=IN+nnz(D);
    pos=pos+nnz(D>0);
    neg=neg+nnz(D<0);
end;
IN;
pos;
neg;
nnz(Slash)-IN;
nnz(Slash>0)-pos;
nnz(Slash<0)-neg;
%tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(Slash,50,0,1,1);toc;%mudulus
%t2=toc;
 
%max(Jf2{2})
 
[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Slash,50,0,1,0);%split
max(Jf1{2})
Q1
 
IN=E<0;
IP=E>0;
n=max(max(X),max(Y));
SSlash=sparse(X(IP),Y(IP),1,n,n);
SSlash=double(SSlash|SSlash');
N=sparse(X(IN),Y(IN),-1,n,n);
clear IN IP 
tic;SSlash(N|N')=-1;toc; %201
[x,y]=find(diag(SSlash));
SSlash=SSlash-sparse(1:n,1:n,diag(SSlash),n,n);
clear N
 
[Jf, Cf, kf, db, Q, acc] = NormalCutS(SSlash,50,2,0);
 
[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(SSlash,50,2,0);
 
[max(Jf1{2}) max(Jf{1}) max(Jf0{1})]
max(Q1)
Q(max(Jf{1})-1)
Q0(max(Jf0{1})-1)
 
 
 
%% Wiki Signed
% I1=[1:7:size(wikisigned,1)];
% I2=I1+1;
% I3=I1+2;% res
% I4=I1+3;% vote
% SRCTAR=wikisigned([I1 I2]);
% % build ID and index
% [ID,~,Idx] = unique(SRCTAR,'stable');
% X=Idx(1:(size(Idx,1)/2));
% Y=Idx(((size(Idx,1)/2)+1):size(Idx,1));
% 
% % extract votes from cell array strings
% CellEvote=wikisigned(I3);
% CellEres=wikisigned(I4);
% Evote=zeros(size(CellEvote,1),1);
% Eres=zeros(size(CellEres,1),1);
% 
% for i=1:size(CellEvote,1)
%     if size(CellEvote{i},2)==5
%         Evote(i)=1;
%     else
%         Evote(i)=-1;
%     end;
%     if size(CellEres{i},2)==5
%         Eres(i)=1;
%     else
%         Eres(i)=-1;
%     end;
% end;
 
% build wiki signed matrix
n=max([X;Y]);
Wiki=sparse(X,Y,E.*(X~=Y),n,n);
[Jf3, kf3, db3, Q3, ang3, acc3,Cf3,v31,d31,v32,d32,I3,cand3] = Augmented_ADJCluster(Wiki,50,1,0,0,1);%AugAdj
tic;[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Wiki,50,0,1,0);toc;%split
t1=toc;
tic;[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v21,d21,v22,d22,I2,cand2] = Generalized_ADJCluster(Wiki,50,0,1,1);toc;%mudulus
t2=toc;
max(Jf1{2})
max(Jf2{2})
 
IN=0;
neg=0;
pos=0;
% detailed information of clustering results
for i=1:max(Jf1{2})
    I=(Jf1{2}==i);
    nnz(I)
    D=Wiki(I,I);
    nnz(D>0)
    nnz(D<0)   
    nnz(D)/(nnz(I)^2)
    IN=IN+nnz(D);
    pos=pos+nnz(D>0);
    neg=neg+nnz(D<0);
end;
IN;
pos;
neg;
nnz(Wiki)-IN;
nnz(Wiki>0)-pos;
nnz(Wiki<0)-neg;
 
 
[Jf1, kf1, db1, Q1, ang1, acc1,Cf1,v11,d11,v12,d12,I1,cand1] = Generalized_ADJCluster(Wiki,50,0,1,0);%split
max(Jf1{2})
Q1
 
IN=E<0;
IP=E>0;
n=max(max(X),max(Y));
SWiki=sparse(X(IP),Y(IP),1,n,n);
SWiki=double(SWiki|SWiki');
N=sparse(X(IN),Y(IN),-1,n,n);
clear IN IP 
tic;SWiki(N|N')=-1;toc; %201
[x,y]=find(diag(SWiki));
SWiki=SWiki-sparse(1:n,1:n,diag(SWiki),n,n);
clear N
 
[Jf, Cf, kf, db, Q, acc] = NormalCutS(SWiki,50,2,0);
 
[Jf0, kf0, db0, Q0, ang0, acc0] = AdjCluster(SWiki,50,2,0);
 
[max(Jf1{2}) max(Jf{1}) max(Jf0{1})]
max(Q1)
Q(max(Jf{1})-1)
Q0(max(Jf0{1})-1)
 


