function [IND, dbt, ang, C] = ADJ(V,c)
% ADJ
% V: input vectors
% c: threshhould of the norm for each vector

[m,n]=size(V);
%get the L2 norm of row vectors
N = sqrt(sum(V.*V,2));
% find the location with norm greater than c
x = find(N>c);
% find the location with norm less than c
y = find(N<=c);

% normalize vectors to the unit n-sphere
for i = 1:n
V(:,i)=V(:,i)./N;
end;

U=V(x,:);

[~,idx] = max(abs(U));

S=zeros(1,n);
for i = 1:n
    S(i)= sign(U(idx(i),i));
end;

% get the initial seed coordinates for K-means
S = diag(S);

% K-means using default distance measure
%  [IND, C] = kmeans(U,n,'Distance','cityblock','start',S,'emptyaction','drop'); 
% SumInnerCD is the vector containing the sum of distances from nodes to the centroid of each cluster
% CDis the matrix containing the distances of each node to each centroid
[IDX, C,~, CD] = kmeans(U,n,'start',S,'emptyaction','drop');

% calculate the Davies-Bouldin-Index
[dbt,ang] = DBI(IDX,C,CD);

IND = zeros(m,1);
IND(x) = IDX;
clear CD

%find the centroid with which the nodes below the threshold have the
%smallest angle, and assign the nodes to the corresponding group
if c>0
    y = find(IND == 0);
    for j = 1:n
        C(j,:) = C(j,:)/norm(C(j,:));
    end;
    for i = 1:size(y,1)
        [~,my] = max(V(y(i),:)*C');
        IND(y(i)) = my;
    end;

end

