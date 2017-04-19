function [Ap,A, partition] = KpartiteGenerator2(k)
% Generate a k partite graph and output adj matrix with partition info
% Input commu size k with k(i) is the number of nodes in this community and
% connection probability matrix p with p(i,j) as the prob between commu i 
% and commu j for j>i and 0 else where
% q please see the paper for definiton
% 

kn = max(size(k));%# of comm
% if kn~=max(size(p))
%     error('Not enough information');
% end;
d=[1.3 1.3];
K= [0,cumsum(k)];
n = K(kn+1);
partition = sparse(n,kn);
for i =1:kn
    i
   
    partition=partition+sparse((K(i)+1):(K(i+1)),i,1,n,kn);
   
end;
partition = logical(partition);

Ap = sparse(n,n);
A=Ap;
% for i = 1:kn
%         warning('iteration for Inner Rand')
%         A(partition(:,i),partition(:,i)) = triu(sprand(k(i),k(i),sum(p(i,:))*0.8),1);
%     for j = (i+1):kn
%         %Temp=RandUS(k(i)+k(j),p(i,j));
%         %Temp(1:k(i),(k(i)+1):(k(i)+k(j)))
%         warning('iteration for Inter Rand')
%         i 
%         j
%         At(partition(:,i),partition(:,j)) = -sprand(k(i),k(j),p(i,j));
%        % At=At+sparse(partition(:,i),partition(:,j),rand(1,k(i)*k(j))<p(i,j),n,n);
%     end;
% end;

for i = 1:kn
        warning('PLUS2 on diagonal')
        i
        tic;
        R = PLUS2(d(1),5500,13);
        %R = PLUS2(d(1),9000,15);
        A(K(i)+1:K(i+1),K(i)+1:K(i+1)) = R(1:k(i),1:k(i));
        toc;
        
        for j = (i+1):kn
             warning('PLUS2 on off diagonal')
            j
            tic;
%            PL=PLUS2(d(2),3300,9);
%            PL=PL(1:k(i),1:k(j));
%            PL=PLUS3(d(2),round(min(k(i),k(j))/250),5,k(i),k(j));
            PL=PLUS3(d(2),2000,10,k(i),k(j));
            Ap((K(i)+1):K(i+1),(K(j)+1):K(j+1))= -PL;
            toc;
            %B((K(j)+1):K(j+1),(K(i)+1):K(i+1))= p(1:k(i),(k(i)+1):(k(i)+k(j)))';
        end;
end;
warning('symmetrization')


Ap = Ap+Ap';

%warning('Eigen decompostion and get q')    
% tic;
% [V,D]=eigs(Ap,kn);
% toc;
%q = zeros(size(p));
% for i = 1:kn
%     for j = (i+1):kn
%         i
%         j
%          warning('iteration for Q')
%         vi = V(partition(:,i),1);
%         vj = V(partition(:,j),1);
%         q(i,j) = D(1,1)*norm(vi)^2*norm(vj)^2/(vi'*Ap(partition(:,i),partition(:,j))*vj)-1;
%        
%     end;
% end;
