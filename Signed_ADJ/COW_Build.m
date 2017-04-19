%Construct COW dataset
%Original Input data
clear
%%
Al = data(:,1:5);
Disp = data(:,6:11);

[n1,n2]=size(Al);
[m1,m2]=size(Disp);


At = sparse(max(Al(:,1)),max(Al(:,2)));
Bt = sparse(max(Disp(:,1)),max(Disp(:,2)));
%% 
for i = 1:n1
    if (Al(i,3)>=1993)
        At(Al(i,1),Al(i,2)) = Al(i,4);
    end;
end;
A = At'*At;
A = A - diag(diag(A));
A = double(logical(A>0));
for i = 1:m1
    if (Disp(i,3)>=1993)
        if (Disp(i,4) == 1)
            Bt(Disp(i,1),Disp(i,2))=Disp(i,5);
        else
            Bt(Disp(i,1),Disp(i,2))=-Disp(i,5);
        end;
    end;        
end;
B = zeros(max(Disp(:,2)),max(Disp(:,2)));
for j = 1:(max(Disp(:,1)))
    if(sum(abs(Bt(j,:)))~=0)
    Bt1 = find(sign(Bt(j,:))==1);
    Bt2 = find(sign(Bt(j,:))==-1);
    B(Bt1,Bt2) = -1;
    end;
end;

B =-double(logical( B + B'));

%%
% B = Bt'*Bt;
% B = B - diag(diag(B));
% B = -double(logical(B<0));
% 
% BD = double(logical(B>0))-double(logical(B<0));

[l1,l2]=size(A);
[l3,l4]=size(B);

if (l3>l1)
    C = B;
    C(1:l1,1:l1) = C(1:l1,1:l1) + A;
else
    C = A;
    C(1:l3,1:l3) = C(1:l3,1:l3) + B;
end;

if (l3>l1)
    C1 = -B;
    C1(1:l1,1:l1) = C1(1:l1,1:l1) + A;
else
    C1 = A;
    C1(1:l3,1:l3) = C1(1:l3,1:l3) - B;
end;
C1 = double(logical(C1));
dif = nnz(C1)-nnz(C);
x = find(sum(abs(C))~=0);
E = C(x,x);
[F,z1]=ExSubG(E,1);

% y1 = find(sum(A)~=0);
% AE = A(y1,y1);
% [d,dm]=max(sum(A));
% AF = ExSubG(AE,dm);

% y2 = find(sum(BD)~=0);
% BE = B(y2,y2);
% [d,dm]=max(sum(A));
% BF = ExSubG(BE,dm);

%%
[xd,yd]=find(triu(C1-abs(C))==1);
[xi,xj]=size(xd);
for i = 1:xi
    if (xd(i)>l1) || (yd(i)>l2) 
        t1 = 0;
    else
        s1 = logical(At(:,xd(i)).*At(:,yd(i)));
        if (sum(s1)==0)
            t1 = 0;
        else
        t1 = sum(s1.*(At(:,xd(i))+At(:,yd(i)))/2)/sum(s1);
        end;
    end;
        s2 = logical(Bt(:,xd(i)).*Bt(:,yd(i))==-1);
        if (sum(s2)==0)
            t2 = 0;
        else
        t2 = sum(s2.*(abs((Bt(:,xd(i)))+abs(Bt(:,yd(i))))/2))/sum(s2);
        end;
    if ((t1<3)&&(t2<3))
        C(xd(i),yd(i))=1;
        C(yd(i),xd(i))=1;
    elseif ((t1>2)&&(t2>2))
        C(xd(i),yd(i))=-1;
        C(yd(i),xd(i))=-1;
    end;
end;

x = find(sum(abs(C))~=0);
CE = C(x,x);
[CF,z2]=ExSubG(CE,1);

