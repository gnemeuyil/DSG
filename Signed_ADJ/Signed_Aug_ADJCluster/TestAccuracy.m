function [ det, mod, acc] = TestAccuracy( G,Jf0,alpha )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[Jfa, kfa, dba, Qa, anga, acca,Cfa] = AdjCluster(G,10,2,0);
%Aug_Adj
[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G,10,2,0,0,alpha);
[Jfs, kfs, dbs, Qs, angs, accs,Cfs] = AdjClusterSVD(G,10,2,0);

det=[max(Jfa{2}) max(Jfaa{2}) max(Jfs{2})]

mod=[max(Qa),max(Qaa),max(Qs)]

Pa=idx2lgc(Jfa{2});
Paa=idx2lgc(Jfaa{2});
Ps=idx2lgc(Jfs{2});
for i=max(Jfa{2}):max(Jf0{2})-1
    Pa=[Pa zeros(size(G,2),1)];
end

for i=max(Jfaa{2}):max(Jf0{2})-1
    Paa=[Paa zeros(size(G,2),1)];
end

for i=max(Jfs{2}):max(Jf0{2})-1
    Ps=[Ps zeros(size(G,2),1)];
end



acc= [PartitionAccuracy( Pa,idx2lgc(Jf0{2}) )
 PartitionAccuracy(Paa,idx2lgc(Jf0{2}) )
PartitionAccuracy(Ps,idx2lgc(Jf0{2}) )]












end

