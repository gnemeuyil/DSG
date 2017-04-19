function [Jf1, kf1, db1, Q1, ang1, acc1,Cf,v1,d1,v2,d2,I,cand] = Generalized_ADJCluster(A,K,c,alpha,MorR,Partition,CL)
% Generalized_ADJCcluster Summary of this function goes here
%   Generalized ADJ cluster methods for directed signed graphs.
%   A:  A square GADJacency matrix of a given graph: directed/undirected
%       signed/unsigned with real edge weights.
%   K:  Desired eigenpairs to probe
%   c:  [0,1] distance threshold for assigning nodes to clusters      
%   alpha:  [0,1] parameter to GADJust the effect of given objective 
%       function. This factor will cause the clustering process to produce 
%       finer or coarser topologies. For example, with smaller alpha 
%       values,clusters could be further divided into smaller clusters. It
%       could help to explore the micro level of graph structures
%   MorR:   0,1 Handle the complex eigenpairs by using Moduli or Splitting
%       method. The moduli mehtod uses the modulus of each entry of the
%       complex vectors multiplied by the sign of its real part. The
%       splitting method splits the complex eigenvectors into real and
%       imaginary parts and map the two value sets into two real vectors
%   CL: Number of clusters is known or expected. Otherwise, ignore it
%   Partition:  True partition results for comparison. Ignore if not known
%   
 
m = size(A,1);
if ~exist('K','var') && ~exist('L','var')&& ~exist('c','var')
    K = min(100, ceil(m/10));
    L = 2;
    c = 0;
end;
    p = 'off';
   % p='off';   
% calculate the K largest positive real eigenvalues with corresponding
% eigenvectors which are all positive
% For signed graphs, use Moduli, Real parts or split Complex eigenpairs 
% into real and imaginary parts by using MorR=1, 2 or 0
tic
if MorR==1
    [V, D1,V2,D2,In ]=ScreenEigpairsSignedC2M(A,K);
elseif MorR==2
    [V, D1,V2,D2,In ]=ScreenEigpairsSignedC2R1(A,K);
else
    [V, D1,V2,D2,In ]=ScreenEigpairsSignedC2R(A,K);
end;
toc
warning('Screening Eigen pairs done!)');
% clustering using K-means
%[Jt, dbt, ang3, Ct] = GADJ(ax,c);
 
% In contains all the location of no-positive/non-negative vectors. Those 
% are part of the vectors we are looking for, Find the candidate vectors 
% by removing all the sets that do not include In.
%
%
% combos=combntns(1:size(V,2),C);
% Lia=ismember(combos,In);
% Insize=size(In,2);
% Loc=sum(Lia,2)==Insize;
% 
% candidates=combos(Loc,:);
% 
% 
%  
% db = zeros(size(candidates,1),3);
% kf = zeros(3,1);
% Q = zeros(size(candidates,1),1);
% ang = zeros(size(candidates,1),3);
% acc = zeros(size(candidates,1),1);
% temp = [10 0];
 
I=1:size(V,2);
CandLoc=ismember(I,In);
NotIn=I(~CandLoc);
Cand=In;
 
 
 
 
 
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
        [Jt, dbt, ang3, Ct] = GADJ(V(:,Cand),c,CL);
    else
        [Jt, dbt, ang3, Ct] = GADJ(V(:,Cand),c);
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
            [Jt, dbt, ang3, Ct] = GADJ(V(:,[Cand NotIn(i-1)]),c,CL);
        else
            [Jt, dbt, ang3, Ct] = GADJ(V(:,[Cand NotIn(i-1)]),c);
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
        [Jt, dbt, ang3, Ct] = GADJ(V(:,[Cand NotIn(i-1)]),c);
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
%     
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
%kf(2)
%size(V)
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
 


