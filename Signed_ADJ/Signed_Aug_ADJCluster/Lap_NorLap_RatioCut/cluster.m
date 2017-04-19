function [ idx,C, sum ] = cluster( L, k )
[V,~]=eigs(L,k);
[n1,n] = size(V);
for i=1:n,
row = V(i,:);
norm_row = norm(row);
for j=1:k,
V(i,j)=V(i,j)/norm_row;
end
end


[~,idx] = max(abs(V));
S=zeros(1,n);
for i = 1:n
    S(i)= sign(V(idx(i),i));
end;

% get the initial seed coordinates for K-means
S = diag(S);
[idx,C, sum]=kmeans(V,k,'start',S,'emptyaction','drop');

%fid = fopen(¡¯partition.txt¡¯,¡¯wt¡¯);
%for i=1:n1,
%fprintf(fid, ¡¯%d\n¡¯, idx(i));
%end
%fclose(fid);
end


