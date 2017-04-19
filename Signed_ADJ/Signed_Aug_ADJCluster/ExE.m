function [Jf,v1,d1,v2,d2,I,cand] = ExE( X,K,c,Shade,pl,C )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[Jf,~, ~, ~,~, ~,~,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(X,K,2,c,C);
if pl==1
%view(biograph(X));
plotClusters(Jf,v,Shade);
end;

end

