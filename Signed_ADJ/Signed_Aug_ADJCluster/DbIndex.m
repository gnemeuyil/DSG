function [dbt1 ang3] = DbIndex(V, Jt1, Ct1, dist_type)
[mC,nC]=size(Ct1);
pt=[];
pho=[];
the=[];

%calculates the distance between and pair of centroids and copy the values
%into the lower half of the matrix to make it symmetric. The trace are 0s
if dist_type==1
    C_dist_M = zeros(mC,mC);
    for s = 1:mC
        for t = s+1:mC
            C_dist_M(s,t)=sum((Ct1(s,:)-Ct1(t,:)).^2);
        end
    end
    C_dist_M = sqrt(C_dist_M+C_dist_M');
else
    C_dist_M = 1-Ct1*Ct1';
end

for s = 1:mC
    for t = [1:(s-1) (s+1):mC]
        u = find(Jt1 == s);
        v = find(Jt1 == t);
        a = mdis(V(u,:),Ct1(s,:),dist_type);
        b = mdis(V(v,:),Ct1(t,:),dist_type);     
        dc = C_dist_M(s,t); 
        % pt calcultes is points belong to different group are closer to
        % their centroids or closer to each other, smaller value indicates
        % points are more closer to thier centroids.
        pt= [pt (a+b)/dc];
        the = [the vectang(Ct1(s,:),Ct1(t,:))];
    end;
    p = max(pt);
    pt = [];
    % records the most probable pair of groups such that the centroids are
    % relatively close but points in each group on average are further away
    % from their centroids for each group
    pho = [pho p];
end;

dbt1 = [mean(pho(~isnan(pho))),max(pho(~isnan(pho))),min(pho(~isnan(pho)))];
ang3 = [mean(the(~isnan(the))),max(the),min(the)];