
[ G_7,~,~] = CombineComponents( G8,G7,0.1,0.1);
[ G_7,~,~] = CombineComponents( G_7,G6,0.1,0.1);
[ G_7,~,~] = CombineComponents( G_7,G5,0.2,0.2);

[ G_7,~,~] = CombineComponents( G_7,G4,0.2,0.2);
[ G_7,~,~] = CombineComponents( G_7,G3,0.2,0.2);
[ G_7,~,~] = CombineComponents( G_7,G1,0.1,0.2);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(G_7,10,2,0,0);
plotClusters( Jf1,v1,[]);
el=adj2edgeL(G_7);
csvwrite('elG7_01_02.csv',el);


[ G_7,~,~] = CombineComponents( G8,G7,0.3,0.3);
[ G_7,~,~] = CombineComponents( G_7,G6,0.5,0.5);
%[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(G_7,10,2,0,0);

[ G_7,~,~] = CombineComponents( G_7,G5,0.2,0.2);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(G_7,10,2,0,0);


[ G_7,~,~] = CombineComponents( G_7,G4,0.2,0.1);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(G_7,10,2,0,0);

[ G_7,~,~] = CombineComponents( G_7,G3,0.2,0.2);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(G_7,10,2,0,0);

[ G_7,~,~] = CombineComponents( G_7,G1,0.1,0.1);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(G_7,10,2,0,0);

[ Jf, Cf, kf, db, Q, acc] = Wcut( G_7,4,2,0,0);
[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v1_star,d1_star,v2_star,d2_star,I2,cand2] = Augmented_ADJCluster(G_star,4,2,0,0);
