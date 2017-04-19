% DBI calculates the Davies-Bouldin-Index
% U: is the coordinates of input vectors
% IDX: contains the cluster indicators of each point
% C: is the coordinates of the centroids
% SumInnerCD: is the sum of distances of nodes to the centroid in a cluster
% CD: the matrix containing the distances of each node to each centroid
%
function [dbt1 ang3] = DBI(IDX,C, CD)

[mC,nC]=size(C);
pt=[];
pho=[];
the=[];
MDIS=[];


%for each cluster calculate the average&max&min distances from the nodes to the centroid 
for s = 1:mC
    s
   u = IDX == s;
   CD_st=CD(u,s);
   mean(CD_st)
   max(CD_st)
   min(CD_st)
   MDIS=[MDIS; mean(CD_st) max(CD_st) min(CD_st)] ;
end;   

 
for s = 1:mC
    for t = [1:(s-1) (s+1):mC]
        dc = sqrt(sum((C(s,:)-C(t,:)).^2)); 
        % pt calcultes is points belong to different group are closer to
        % their centroids or closer to each other, smaller value indicates
        % points are more closer to thier centroids.
        pt= [pt (MDIS(s,:)+MDIS(t,:))/dc];
        the = [the vectang(C(s,:),C(t,:))];
    end;
    p = max(pt);
    pt = [];
    % records the most probable pair of groups such that the centroids are
    % relatively close but points in each group on average are further away
    % from their centroids for each group
    pho = [pho p];
end;

dbt1 = [mean(pho(~isnan(pho))),max(pho(~isnan(pho))),min(pho(~isnan(pho)))];
ang3 = [min(the),mean(the(~isnan(the))),max(the)];
clear u the pho MDIS pt s t
