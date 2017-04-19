function [ Jf1,Q1,v1,d1,I,cand ] = Sym( G1,G2,p1,p2,K )
%SYM Summary of this function goes here
%   Detailed explanation goes here

G=CombineComponents(G1, G2,p1,p2);


SG=G|G';
SG=double(SG);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(SG,K,2,0,0);
plotClusters( Jf1,v1 ,[]);
end

