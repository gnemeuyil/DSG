function [Jt1, dbt1, ang3, Ct] = WAdjCl(V,L,c)
[m,n]=size(V);
N = sqrt(sum(V.*V,2));
x = N>c;

U=V;
[~,idx] = max(abs(U));

S=zeros(1,n);
for i = 1:n
    S(i)= sign(U(idx(i),i));
end;
S = diag(S);
% S = diag(S);
% S = [S;zeros(1,n)];
% S(n+1,1)=-S(1,1);
if L==1
    [Jt, Ct] = kmeans(U,n,'Distance','cityblock','start',S,'emptyaction','drop'); 
elseif L==2
    [Jt, Ct] = kmeans(U,n,'start',S,'emptyaction','drop'); 
    %[Jt Ct] = kmeans(U,n+1);
end;

[dbt1,ang3] = DbIdx(U,Jt,Ct,1);

Jt1 = zeros(m,1);
Jt1(x) = Jt;

if c>0
    y = find(Jt1 == 0);
    for j = 1:n
        Ct(j,:) = Ct(j,:)/norm(Ct(j,:));
    end;
    for i = 1:size(y,1)
        [~,my] = max(V(y(i),:)*Ct');
        Jt1(y(i)) = my;
    end;
end;

