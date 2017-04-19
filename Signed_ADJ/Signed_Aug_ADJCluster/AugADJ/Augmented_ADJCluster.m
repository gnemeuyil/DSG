function [Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Augmented_ADJCluster(A,K,L,c,C,alpha,Partition)
%AUGMENTED_ADJCLUSTER Summary of this function goes here
%   Detailed explanation goes here

m = size(A,1);
if ~exist('K','var') && ~exist('L','var')&& ~exist('c','var')
    K = min(100, ceil(m/10));
    L = 2;
    c = 0;
end;
    p = 'off';
   % p='off';
   
   
   tic;



% calculate the K largest positive real eigenvalues with corresponding
% eigenvectors which are all positive
%[ax d1]=FindEigpairs(A,K);
tic
[V, D1,V2,D2,In ]=ScreenEigpairs(A,K);
toc
warning('Screening Eigen pairs done!)');





I=1:size(V,2);
CandLoc=ismember(I,In);
NotIn=I(~CandLoc);
Cand=In;

% db = zeros(size(NotIn,2)+1,3);
% kf = zeros(3,1);
% Q = zeros(size(NotIn,2)+1,1);
% ang = zeros(size(NotIn,2)+1,3);
% acc = zeros(size(NotIn,2)+1,1);
% temp = [10 0];
% 
% 
% % calculate Modularity and DBI for initial candidate set 
%    tic
%   
%     [Jt, dbt, ang3, Ct] = ADJ(V(:,Cand),c);
%     
%     db(1,:) = dbt;
%     ang(1,:) = ang3;
%     if db(1,1)<=temp(1)
%         temp(1) = db(1,1);
%         Jf{1} = Jt;
%         kf(1) = size(Cand,2);
%         Cf{1} = Ct;
%     elseif temp(1) == db(1,1)
%         warning('Multiple Solution!)');
%     end;
%     
%     
%     
%     P = idx2lgc(Jt);
%     Q(1) = SignQfunction(A,P);
% 
%     if Q(1)>=temp(2)
%         temp(2) = Q(1);
%         Jf{2} = Jt;
%         kf(2) = size(Cand,2);
%         Cf{2} = Ct;
%     elseif temp(2) == Q(1)
%        % warning('Multiple Solution!)');
%     end;
%     
%     if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
%         acc(1) = PartitionAccuracy(P, Partition);
%     end
%     
%  
%      warning('Initialize Modularity and DBI with candidate set)');
%    toc
% 
% 
% 
% 
% 
% 
% 
% 
% 
% CandDBI=Cand;
% 
% 
% 
% % combinatorial fitting based on modularity measures
% for i = 2:size(NotIn,2)+1
%     tic
%     i
%     [Jt, dbt, ang3, Ct] = ADJ(V(:,[Cand NotIn(i-1)]),c);
%     
%     db(i,:) = dbt;
%     ang(i,:) = ang3;
%     if temp(1)>db(i,1)
%         temp(1) = db(i,1);
%         Jf{1} = Jt;
%         CandDBI=[CandDBI NotIn(i-1)];
%         kf(1) = size(CandDBI,2);
%         Cf{1} = Ct;
%     elseif temp(1) == db(i,1)
%         warning('Multiple Solution!)');
%     end;
%     
%     
%     
%     P = idx2lgc(Jt);
%     Q(i) = SignQfunction(A,P)
% 
%     if Q(i)>temp(2)*alpha
%         
%         temp(2) = Q(i)
%         Jf{2} = Jt;
%        Cand=[Cand NotIn(i-1)]
%         kf(2) = size(Cand,2);
%         Cf{2} = Ct;
%         
%     elseif temp(2) == Q(i)
%         warning('Multiple Solution!)');
%     end;
%     
%     if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
%         acc(i) = PartitionAccuracy(P, Partition);
%     end
%     
%  
%      warning('one Loop done!)');
%    toc
%   
% end;
% 
% 
% toc;
% t=toc;
% calculate Modularity and DBI for initial candidate set if exist
if ~isempty(Cand)
    warning('Cand not empty')
    db = zeros(size(NotIn,2)+1,3);
    kf = zeros(3,1);
    Q = zeros(size(NotIn,2)+1,1);
    ang = zeros(size(NotIn,2)+1,3);
    acc = zeros(size(NotIn,2)+1,1);
    temp = [inf 0];
    if exist('CL','var')&&CL<=nnz(Cand)
        [Jt, dbt, ang3, Ct] = ADJ(V(:,Cand),c,CL);
    else
        [Jt, dbt, ang3, Ct] = ADJ(V(:,Cand),c);
    end;
    db(1,:) = dbt;
    ang(1,:) = ang3;
    if db(1,1)<=temp(1)
        temp(1) = db(1,1);
        Jf{1} = Jt;
        kf(1) = size(Cand,2);
        Cf{1} = Ct;
    elseif temp(1) == db(1,1)
        warning('Multiple Solution!)');
    end;
    
    
    
    P = idx2lgc(Jt);
    Q(1) = SignQfunction(A,P);
 
    if abs(Q(1))>=temp(2)%modularity may be negative for signed graphs
        temp(2) = Q(1);
        Jf{2} = Jt;
        kf(2) = size(Cand,2);
        Cf{2} = Ct;
    elseif temp(2) == Q(1)
       % warning('Multiple Solution!)');
    end;
    
    if exist('Partition','var') && size(P,2)-size(Partition,2)>=0;
        PartitionE=zeros(size(P));
        PartitionE(:,1:size(Partition,2))=Partition;
        acc(1) = PartitionAccuracy(P, PartitionE);
    end
 
    warning('Initialize Modularity and DBI with candidate set)');
    CandDBI=Cand;
    % combinatorial fitting based on modularity measures
    for i = 2:size(NotIn,2)+1
       
        i
        if exist('CL','var')&&CL<=nnz(Cand)
            [Jt, dbt, ang3, Ct] = ADJ(V(:,[Cand NotIn(i-1)]),c,CL);
        else
            [Jt, dbt, ang3, Ct] = ADJ(V(:,[Cand NotIn(i-1)]),c);
        end;
        
 
        db(i,:) = dbt;
        ang(i,:) = ang3;
        if temp(1)>db(i,1)
            temp(1) = db(i,1);
            Jf{1} = Jt;
            CandDBI=[CandDBI NotIn(i-1)];
            kf(1) = size(CandDBI,2);
            Cf{1} = Ct;
        elseif temp(1) == db(i,1)
            warning('Multiple Solution!)');
        end;
 
 
 
        P = idx2lgc(Jt);
        Q(i) = SignQfunction(A,P);
 
        if abs(Q(i))>abs(temp(2))*alpha%modularity may be negative for signed graphs
 
            temp(2) = Q(i);
            Jf{2} = Jt;
            Cand=[Cand NotIn(i-1)];
            kf(2) = size(Cand,2);
            Cf{2} = Ct;
 
        elseif temp(2) == Q(i)
            warning('Multiple Solution!)');
        end;
        %sum((size(P)-size(Partition)).^2)
        if exist('Partition','var') && size(P,2)-size(Partition,2)>=0;
        PartitionE=zeros(size(P));
        PartitionE(:,1:size(Partition,2))=Partition;
        acc(i-1) = PartitionAccuracy(P, PartitionE);
        end
 
    end;   
else
    % if initial candidate set is empty, iterate through all eigen pairs to
    % find an optimal solution
    warning('Cand empty')
    Cand=[Cand NotIn(1)];
    NotIn=NotIn(2:size(NotIn,2));
    CandDBI=Cand;
    
    db = zeros(size(NotIn,2)+1,3);
    kf = zeros(3,1);
    Q = zeros(size(NotIn,2)+1,1);
    ang = zeros(size(NotIn,2)+1,3);
    acc = zeros(size(NotIn,2)+1,1);
    temp = [inf 0];
    
    
      % combinatorial fitting based on modularity measures
    for i = 2:size(NotIn,2)+1
       i
        [Jt, dbt, ang3, Ct] = ADJ(V(:,[Cand NotIn(i-1)]),c);
        max(Jt)
 
        db(i,:) = dbt;
        ang(i,:) = ang3;
        if temp(1)>db(i,1)
            temp(1) = db(i,1);
            Jf{1} = Jt;
            CandDBI=[CandDBI NotIn(i-1)];
            kf(1) = size(CandDBI,2);
            Cf{1} = Ct;
        elseif temp(1) == db(i,1)
            warning('Multiple Solution!)');
        end;
 
 
 
        P = idx2lgc(Jt);
        Q(i) = SignQfunction(A,P);
        if abs(Q(i))>abs(temp(2))*alpha%modularity may be negative for signed graphs
 
            temp(2) = Q(i);
            Jf{2} = Jt;
            Cand=[Cand NotIn(i-1)];
            kf(2) = size(Cand,2);
            Cf{2} = Ct;
 
        elseif temp(2) == Q(i)
            warning('Multiple Solution!)');
        end;
        size(P)
        %size(Partition)
        
        if exist('Partition','var') && size(P,2)-size(Partition,2)>=0;
        PartitionE=zeros(size(P));
        PartitionE(:,1:size(Partition,2))=Partition;
        acc(i-1) = PartitionAccuracy(P, PartitionE);
        end
      
    end;  
    
end





% if strcmp(p,'on')
%     figure;plot(1:size(NotIn,2)+1,db(:,1)','-o','MarkerSize',6,'MarkerFaceColor','b');
%     xlabel('k','FontSize',14);
%     ylabel('Davies-Buldin Index','FontSize',14);
%  
%     box on;
%     axis square;
% 
%     
%     
%     
%     figure;plot(1:size(NotIn,2)+1,abs(Q),'-o','MarkerSize',6,'MarkerFaceColor','b');
%     xlabel('k','FontSize',14);
%     ylabel('Modularity','FontSize',14);
% 
%     box on;
%     axis square;
% 
% 
% end;

% v1 are eigenvectors by modularity v2 are eigenvectors of real positive
% entries
kf(2)
size(V)
v1=V(:,Cand);
d1=D1(Cand);
v2=V2;
d2=D2;
Jf1 = Jf;
kf1 = kf;
db1 = db;
Q1 = Q;
ang1 = ang;
acc1 = acc;
I=In;    
cand=Cand;   
    
    
    
    
end

