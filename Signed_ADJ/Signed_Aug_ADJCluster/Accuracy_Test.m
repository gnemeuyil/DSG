%% 3 communities disconnected

[ G_123_0_0,~,~] = CombineComponents( G3,G2,0,0);
[ G_123_0_0,~,~] = CombineComponents( G_123_0_0,G1,0,0);


[ G_123_04_05,~,~] = CombineComponents( G3,G2,0.5,0.5);
[ G_123_04_05,~,~] = CombineComponents( G_123_04_05,G1,0.1,0.1);

el=adj2edgeL(G_123_04_05);
csvwrite('G_3_2_05_1_01.csv',el);
%Weighted cut
[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_123_04_05,10,2,0,0);
%Normalized cut
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_123_04_05,10,2,0);
%ADJ
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_123_04_05,10,2,0);
%Aug_Adj
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_123_04_05,10,2,0,0);


acc(1) = PartitionAccuracy([zeros(148,1) idx2lgc(Jf{2})],idx2lgc(Jf0{2}) );
acc(2) = PartitionAccuracy([zeros(148,1) idx2lgc(Jfn{2})],idx2lgc(Jf0{2}) );
acc(3) = PartitionAccuracy([zeros(148,1) idx2lgc(Jfa{2})],idx2lgc(Jf0{2}) );
acc(4) = PartitionAccuracy(idx2lgc(Jfaa{2}),idx2lgc(Jf0{2}) );

%% 7 communities disconnected

[ G_7_0_0_18,~,~] = CombineComponents( G8,G7,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G6,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G5,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G4,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G3,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G1,0,0);

%Weighted cut
[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_0_0_18,10,2,0,0);
%Normalized cut
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_0_0_18,10,2,0);
%ADJ
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_0_0_18,10,2,0);
%Aug_Adj
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_0_0_18,10,2,0,0);

%% 7 communities balanced 1
[ G_7_01,~,~] = CombineComponents( G8,G7,0.1,0.1);
[ G_7_01,~,~] = CombineComponents( G_7_01,G6,0.1,0.1);
[ G_7_01,~,~] = CombineComponents( G_7_01,G5,0.1,0.1);
[ G_7_01,~,~] = CombineComponents( G_7_01,G4,0.1,0.1);
[ G_7_01,~,~] = CombineComponents( G_7_01,G3,0.1,0.1);
[ G_7_01,~,~] = CombineComponents( G_7_01,G1,0.1,0.1);

[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_01,10,2,0,0);
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_01,10,2,0);
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_01,10,2,0);
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_01,10,2,0,0)
[Jfs, kfs, dbs, Qs, angs, accs,Cfs] = AdjClusterSVD(G_7_01,10,2,0);
max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})
max(Jfs{2})
acc(1) = PartitionAccuracy(idx2lgc(Jf{2}),idx2lgc(Jf0{2}) )
acc(2) = PartitionAccuracy(idx2lgc(Jfn{2}),idx2lgc(Jf0{2}) )
acc(3) = PartitionAccuracy(idx2lgc(Jfa{2}),idx2lgc(Jf0{2}) )
acc(4) = PartitionAccuracy(idx2lgc(Jfaa{2}),idx2lgc(Jf0{2}) )

%% 7 communities balanced 2 formerly G7_05_6_01_1_01 **
[ G_7_025_02_6_01_01_1_03,~,~] = CombineComponents( G8,G7,0.25,0.2);
[ G_7_025_02_6_01_01_1_03,~,~] = CombineComponents( G_7_025_02_6_01_01_1_03,G6,0.1,0.1);
[ G_7_025_02_6_01_01_1_03,~,~] = CombineComponents( G_7_025_02_6_01_01_1_03,G5,0,0);
[ G_7_025_02_6_01_01_1_03,~,~] = CombineComponents( G_7_025_02_6_01_01_1_03,G4,0,0);
[ G_7_025_02_6_01_01_1_03,~,~] = CombineComponents( G_7_025_02_6_01_01_1_03,G3,0,0);
[ G_7_025_02_6_01_01_1_03,~,~] = CombineComponents( G_7_025_02_6_01_01_1_03,G1,0,0.3);

[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_025_02_6_01_01_1_03,10,2,0,0);
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_025_02_6_01_01_1_03,10,2,0);
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_025_02_6_01_01_1_03,10,2,0);
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_025_02_6_01_01_1_03,10,2,0,0)

max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})

max(Q)
max(Qn)
max(Qa)
max(Qaa)

acc(1) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jf{2})],idx2lgc(Jf0{2}) );
acc(2) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfn{2})],idx2lgc(Jf0{2}) );
acc(3) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfa{2})],idx2lgc(Jf0{2}) );
acc(4) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfaa{2})],idx2lgc(Jf0{2}) );

%% 7 communities unbalanced sparse 1 **
[ G_7_01_02,~,~] = CombineComponents( G8,G7,0.1,0.2);
[ G_7_01_02,~,~] = CombineComponents( G_7_01_02,G6,0.1,0.2);
[ G_7_01_02,~,~] = CombineComponents( G_7_01_02,G5,0.2,0.1);
[ G_7_01_02,~,~] = CombineComponents( G_7_01_02,G4,0.2,0.1);
[ G_7_01_02,~,~] = CombineComponents( G_7_01_02,G3,0.1,0.2);
[ G_7_01_02,~,~] = CombineComponents( G_7_01_02,G1,0.2,0.1);

[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_01_02,10,2,0,0);
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_01_02,10,2,0);
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_01_02,10,2,0);
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_01_02,10,2,0,0)

max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})

acc(1) = PartitionAccuracy(idx2lgc(Jf{2}),idx2lgc(Jf0{2}) )
acc(2) = PartitionAccuracy(idx2lgc(Jfn{2}),idx2lgc(Jf0{2}) )
acc(3) = PartitionAccuracy(idx2lgc(Jfa{2}),idx2lgc(Jf0{2}) )
acc(4) = PartitionAccuracy(idx2lgc(Jfaa{2}),idx2lgc(Jf0{2}) )

%% 7 com unbalanced 2 **
[ G_7_87_01_05_6_03_01_1_05,~,~] = CombineComponents( G8,G7,0.1,0.5);
[ G_7_87_01_05_6_03_01_1_05,~,~] = CombineComponents( G_7_87_01_05_6_03_01_1_05,G6,0.3,0.1);
[ G_7_87_01_05_6_03_01_1_05,~,~] = CombineComponents( G_7_87_01_05_6_03_01_1_05,G5,0,0);
[ G_7_87_01_05_6_03_01_1_05,~,~] = CombineComponents( G_7_87_01_05_6_03_01_1_05,G4,0,0);
[ G_7_87_01_05_6_03_01_1_05,~,~] = CombineComponents( G_7_87_01_05_6_03_01_1_05,G3,0,0);
[ G_7_87_01_05_6_03_01_1_05,~,~] = CombineComponents( G_7_87_01_05_6_03_01_1_05,G1,0,0.5);

[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_87_01_05_6_03_01_1_05,10,2,0,0);
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_87_01_05_6_03_01_1_05,10,2,0);
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_87_01_05_6_03_01_1_05,10,2,0);
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_87_01_05_6_03_01_1_05,10,2,0,0);

max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})

acc(1) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jf{2})],idx2lgc(Jf0{2}) )
acc(2) = PartitionAccuracy( idx2lgc(Jfn{2}),idx2lgc(Jf0{2}) )
acc(3) = PartitionAccuracy( idx2lgc(Jfa{2}),idx2lgc(Jf0{2}) )
acc(4) = PartitionAccuracy( idx2lgc(Jfaa{2}),idx2lgc(Jf0{2}) )

%% Dense 
[ G_7_025_01,~,~] = CombineComponents( G8,G7,0.25,0.25);
[  G_7_025_01,~,~] = CombineComponents( G_7_025_01,G6,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G5,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G4,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G3,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G1,0.1,0.1);

[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_025_01,25,2,0,0);
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_025_01,25,2,0);
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_025_01,25,2,0);
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_025_01,25,2,0,0);

max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})


max(Q)
max(Qn)
max(Qa)
max(Qaa)

acc(1) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jf{2})],idx2lgc(Jf0{2}) )
acc(2) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfn{2})],idx2lgc(Jf0{2}) )
acc(3) = PartitionAccuracy([zeros(1238,1) zeros(1238,1) zeros(1238,1) idx2lgc(Jfa{2})], idx2lgc(Jf0{2}) )
acc(4) = PartitionAccuracy( idx2lgc(Jfaa{2}),idx2lgc(Jf0{2}) )

%% 18 nodes
[ G_7_0_0_18,~,~] = CombineComponents( G8,G7,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G6,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G5,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G4,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G3,0,0);
[ G_7_0_0_18,~,~] = CombineComponents( G_7_0_0_18,G,0,0);
K=10;
[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_0_0_18,K,2,0,0);
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_0_0_18,K,2,0);
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_0_0_18,K,2,0);
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_0_0_18,K,2,0,0);
[Jfs, kfs, dbs, Qs, angs, accs,Cfs] = AdjClusterSVD(G_7_0_0_18,K,2,0);
max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})
max(Jfs{2})

acc(1) = PartitionAccuracy(idx2lgc(Jf{2}),idx2lgc(Jf0{2}) )
acc(2) = PartitionAccuracy(idx2lgc(Jfn{2}),idx2lgc(Jf0{2}) )
acc(3) = PartitionAccuracy([zeros(1228,1) idx2lgc(Jfa{2})], idx2lgc(Jf0{2}) )
acc(4) = PartitionAccuracy( idx2lgc(Jfaa{2}),idx2lgc(Jf0{2}) )
acc(5) = PartitionAccuracy( [zeros(1228,1) idx2lgc(Jfs{2})],idx2lgc(Jf0{2}) )


%% 7 communities disconnected

[ G_7_0_0,~,~] = CombineComponents( G8,G7,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G6,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G5,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G4,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G3,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G1,0,0);

%Weighted cut
[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_0_0,10,2,0,0);
%Normalized cut
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_0_0,10,2,0);
%ADJ
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_0_0,10,2,0);
%Aug_Adj
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_0_0,10,2,0,0);
[Jfs, kfs, dbs, Qs, angs, accs,Cfs] = AdjClusterSVD(G_7_0_0,K,2,0);

max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})
max(Jfs{2})
acc(3) = PartitionAccuracy( [zeros(1238,1) zeros(1238,1) zeros(1238,1) idx2lgc(Jfa{2})],idx2lgc(Jf0{2}) )

acc(5) = PartitionAccuracy( [zeros(1238,1) zeros(1238,1) zeros(1238,1) idx2lgc(Jfs{2})],idx2lgc(Jf0{2}) )

%% 7 communities 7-8 0.5 others 0.1

[Jf0, kf0, db0, Q0, ang0, acc0,Cf0] = AdjCluster(G_7_0_0_28,10,2,0);

[ G_7_025_01,~,~] = CombineComponents( G8,G7,0.5,0.5);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G6,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G5,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G4,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G3,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G1,0.1,0.1);

%Weighted cut
[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_025_01,10,2,0,0);
%Normalized cut
[Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_025_01,10,2,0);
%ADJ
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_025_01,10,2,0);
%Aug_Adj
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_025_01,10,2,0,0,0.9);
[Jfs, kfs, dbs, Qs, angs, accs,Cfs] = AdjClusterSVD(G_7_025_01,10,2,0);

max(Jf{2})
max(Jfn{2})
max(Jfa{2})
max(Jfaa{2})
max(Jfs{2})
acc(1) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jf{2})],idx2lgc(Jf0{2}) )
acc(2) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfn{2})],idx2lgc(Jf0{2}) )
acc(3) = PartitionAccuracy( [zeros(1238,1) zeros(1238,1) idx2lgc(Jfa{2})],idx2lgc(Jf0{2}) )
acc(4) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfaa{2})],idx2lgc(Jf0{2}) )
acc(5) = PartitionAccuracy( [zeros(1238,1) zeros(1238,1) idx2lgc(Jfs{2})],idx2lgc(Jf0{2}) )

%% 7 communities 7-8 0.25 others 0.1

[Jf0, kf0, db0, Q0, ang0, acc0,Cf0,v10,d10,v20,d20,I20,cand20] = Augmented_ADJCluster(G_7_0_0_28,10,2,0,0,1);

[ G_7_025_01,~,~] = CombineComponents( G8,G7,0.25,0.25);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G6,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G5,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G4,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G3,0.1,0.1);
[ G_7_025_01,~,~] = CombineComponents( G_7_025_01,G1,0.1,0.1);

[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_025_01,10,2,0,0,0.9);
[det,mod,acc]=TestAccuracy(G_7_025_01,Jf0,0.9);

% %Weighted cut
% [ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_7_05_01,10,2,0,0);
% %Normalized cut
% [Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(G_7_05_01,10,2,0);
% %ADJ
% [Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G_7_05_01,10,2,0);
% %Aug_Adj
% [Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_05_01,10,2,0,0,0.9);
% [Jfs, kfs, dbs, Qs, angs, accs,Cfs] = AdjClusterSVD(G_7_05_01,10,2,0);
% 
% max(Jf{2})
% max(Jfn{2})
% max(Jfa{2})
% max(Jfaa{2})
% max(Jfs{2})
% acc(1) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jf{2})],idx2lgc(Jf0{2}) )
% acc(2) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfn{2})],idx2lgc(Jf0{2}) )
% acc(3) = PartitionAccuracy( [zeros(1238,1) zeros(1238,1) idx2lgc(Jfa{2})],idx2lgc(Jf0{2}) )
% acc(4) = PartitionAccuracy([zeros(1238,1) idx2lgc(Jfaa{2})],idx2lgc(Jf0{2}) )
% acc(5) = PartitionAccuracy( [zeros(1238,1) zeros(1238,1) idx2lgc(Jfs{2})],idx2lgc(Jf0{2}) )

%% test

[ G_7_0_0,~,~] = CombineComponents( G8,G7,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G6,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G5,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G4,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G3,0,0);
[ G_7_0_0,~,~] = CombineComponents( G_7_0_0,G1,0,0);


Jf0 = Augmented_ADJCluster(G_7_0_0,10,2,0,0,1);
Jf28 = Augmented_ADJCluster(G_7_0_0_28,10,2,0,0,1);

temp1=0;
temp2=0;
for i=1:5
    
    Jf1=Augmented_ADJCluster(G_7_0_0,10,2,0,0,1);
    
    temp1=temp1+sum(Jf1{2}-Jf0{2})
    
    temp2=temp2+sum(Jf1{2}-Jf28{2})
    
end


