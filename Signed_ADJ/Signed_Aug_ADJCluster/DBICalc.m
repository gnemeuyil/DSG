function [ db,ang, Jf,kf, Cf ] = DBICalc(ax,K, c,db,ang,acc, Jf,kf,Cf, temp,Q)
%DBICALC Summary of this function goes here
%   Detailed explanation goes here
for i = 2:min(d,K)
   
    [Jt, dbt, ang3, Ct] = AdjCl(ax(:,1:i),L,c);
    
    db(i-1,:) = dbt;
    ang(i-1,:) = ang3;
    if db(i-1,1)<temp(1)
        temp(1) = db(i-1,1);
        Jf{1} = Jt;
        kf(1) = i;
        Cf{1} = Ct;
    elseif temp(1) == db(i-1,1)
        WARNING('Multiple Solution!)');
    end;
    
    
    
    P = idx2lgc(Jt);
    Q(i-1) = SignQfunction(A,P);
%     if strcmp(s, 'lm') && min(min(A))==0
%         Q(i-1)=-Q(i-1);
%     end;
    if Q(i-1)>temp(2)
        temp(2) = Q(i-1);
        Jf{2} = Jt;
        kf(2) = i;
        Cf{2} = Ct;
    elseif temp(2) == Q(i-1)
        WARNING('Multiple Solution!)');
    end;
    
    if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
        acc(i-1) = PartitionAccuracy(P, Partition);
    end
end;
toc;
t=toc;




end

