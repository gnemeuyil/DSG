function [IND, dbt, ang, C] = GADJ(V,c, CL)
% ADJ
% V: input vectors
% c: threshhould of the norm for each vector
% CL: real cluster number

[m,n]=size(V);
loc=[];
for i=1:n
    loc(i)=isreal(V(:,i));%find real vector locations
end;

%for each complex vector split it into to vectors containing real and
%imaginary parts respectively
loc=logical(loc);
cloc=~loc;
U=[V(:,loc) real(V(:,cloc)) imag(V(:,cloc))];

%get the L2 norm of row vectors
N = sqrt(sum(U.*U,2));
% find the location with norm greater than c
x = find(N>c);
% find the location with norm less than c
%y = find(N<=c);
[nu,mu]=size(U);
% normalize vectors to the unit n-sphere
for i = 1:mu
    U(:,i)=U(:,i)./N;
end;

U=U(x,:);

[~,idx] = max(abs(U));



S=zeros(1,n);
for i = 1:n
    S(i)= sign(U(idx(i),i));
end;

% get the initial seed coordinates for K-means
S = sparse(1:n,1:n,S,n,mu);
S=full(S);
% K-means using default distance measure
%  [IND, C] = kmeans(U,n,'Distance','cityblock','start',S,'emptyaction','drop'); 
% SumInnerCD is the vector containing the sum of distances from nodes to the centroid of each cluster
% CDis the matrix containing the distances of each node to each centroid
if exist('CL','var')
    [IDX, C,~, CD] = kmeans(U,CL,'emptyaction','drop'); 
else
    [IDX, C,~, CD] = kmeans(U,n,'start',S,'emptyaction','drop'); 
end;
% calculate the Davies-Bouldin-Index

% remove NAN centroids and distances 
C=C(~any(isnan(C),2),:);
CD=CD(:,~any(isnan(CD),1));
UIDX=unique(IDX);
% shrink indices accordingly
for i=1:size(UIDX,1)
    IDX(IDX==UIDX(i))=i;
end;

[dbt,ang] = DBI(IDX,C,CD);

IND = zeros(m,1);
IND(x) = IDX;


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

