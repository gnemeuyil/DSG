function [ D ] = Distmatrix( A)
%

[n,~]=size(A);
D=sparse(zeros(n,n));
for i=1:n
    tic
for j=1:n
    D(i,j)=dijkstra(A,i,j);
end
i
toc
end



end

