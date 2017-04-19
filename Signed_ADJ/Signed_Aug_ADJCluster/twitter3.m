files1 = dir('Edge_*.txt');
filenames1={files1.name};
[SF1 C1]=sort_nat(filenames1);

Edge=[];
for i=1:length(files1)
    Edge=[Edge; load(char(SF1(i)))];
    SF1(i)
end

files = dir('Index_*.txt');
filenames={files.name};
[SF C]=sort_nat(filenames);

IN=[];
for i=1:length(files)
    IN=[IN; load(char(SF(i)))];
    SF(i)
end;
Index=IN(:,1);


size(Edge)
m=max(max(Edge))
val=ones(size(Edge,1),1);
A=sparse(Edge(:,1),Edge(:,2),val,m,m);
clear Edge SF1 ans files1 filenames1
clear IN files

A2=double(logical(A*A'));
A1=A+A';
[AN1,I1,n1,e1]=precheck(A1);
dlmwrite('I1.csv', find(I1==1)', 'delimiter', ',', 'precision', 10);
dlmwrite('n1.csv', n1, 'delimiter', ',', 'precision', 10);
dlmwrite('e1.csv', e1, 'delimiter', ',', 'precision', 10);
%dlmwrite('AN1.csv', AN1, 'delimiter', ',', 'precision', 10);


AN1=AN1-diag(diag(AN1));

k=20;
[V,D]=eigs(AN1,k);
[Jf1, kf1, db1, Q1, ang1, acc1,Cf1] = AdjCluster(AN1,k,2,0);
%dlmwrite('Jf1.csv', find(Jf1==1)', 'delimiter', ',', 'precision', 10);
dlmwrite('kf1.csv', kf1, 'delimiter', ',', 'precision', 10);
dlmwrite('db1.csv', db1, 'delimiter', ',', 'precision', 10);
dlmwrite('Q1.csv', Q1, 'delimiter', ',', 'precision', 10);
dlmwrite('ang1.csv', ang1, 'delimiter', ',', 'precision', 10);
dlmwrite('acc1.csv', acc1, 'delimiter', ',', 'precision', 10);
%dlmwrite('Cf1.csv',Cf1, 'delimiter', ',', 'precision', 10);
dlmwrite('eigVal.csv',D, 'delimiter', ',', 'precision', 10);
dlmwrite('eigVecl.csv',V, 'delimiter', ',', 'precision', 10);
k=18;
plotopt(AN1,k);


[Jt1 dbt1 ang3 Ct1] = AdjCl(V(:,1:k),2,0);

dlmwrite('dbt1.csv', dbt1, 'delimiter', ',', 'precision', 10);
dlmwrite('ang3.csv', ang3, 'delimiter', ',', 'precision', 10);
dlmwrite('Ct1.csv', Ct1, 'delimiter', ',', 'precision', 10);
dlmwrite('V1_k.csv', V(:,1:k), 'delimiter', ',', 'precision', 10);
n = size(V,1);
[mC,nC]=size(Ct1);
% N is the distance of a node to the origin,bigger number indicates that node is more active/significant 
N = sqrt(sum(V(:,1:k).*V(:,1:k),2));
dlmwrite('N.csv', N, 'delimiter', ',', 'precision', 10);
dlmwrite('Jt1.csv', Jt1, 'delimiter', ',', 'precision', 10);


figure;bar(N,'DisplayName','N');



%% page rank
At=logical(A');
RS=sum(At,2);
p1=find(RS>0);
RS1=RS(p1);
RS2=1./RS1;
W=sparse(p1,p1,RS2,m,m);
WA=W*At;
sum(WA(1:10,:),2)


PR=ones(size(WA,1),1);

PR=PR/m;

Pos_0=(RS==0)/m;
sum(Pos_0)

 for i=1:30
PR=0.15/m+0.85*WA'*PR +0.85*Pos_0'*PR;
i
end;
dlmwrite('PR.csv', PR, 'delimiter', ',', 'precision', 10);
SPR=sort(PR,'descend');
dlmwrite('SPR.csv', SPR, 'delimiter', ',', 'precision', 10);
figure;bar(PR,'DisplayName','PR');

%% results
topIdx=find(PR>=SPR(i));
% or top5=find(ismember(PR,SPR(1:5)));
ID=Index(topIdx);

dlmwrite('Top_PR_ID.csv', [ID PR(topIdx)], 'delimiter', ',', 'precision', 10);
[topN topPR]=Group_Top(N, PR, Jt1,5,Index);
dlmwrite('top_each_group.csv', [topN; topPR], 'delimiter', ',', 'precision', 10);



corporations=[19336500;361259570;81192849;19617284;52522194;127975950;54060157;
    23261531;30981244;20862686;40879784;95042315;15855468;11399812;17325073;26858880;
    70979930;50687788;204856343;101650367;207507805;204881628;18735040;111594658;
    87777436;1178011;23002858;17137891;144841605;89084561;80374332;79320096;426159377;
    16676528;281285283;18083344;20683724;1224620941;20187750;23450421;487624974;34685994];
[corporations_found corporationsIdx group SpectralScore SpectralScoreRank PRScore PRScoreRank position below positionPR belowPR total indicator]=NodeRank(corporations, Index, Jt1,N,PR);
results=[corporations_found corporationsIdx group SpectralScore SpectralScoreRank PRScore PRScoreRank position below positionPR belowPR total];
dlmwrite('results.csv', results, 'delimiter', ',', 'precision', 10);
dlmwrite('corporations_indicator.csv', [corporations indicator], 'delimiter', ',', 'precision', 10);


%% compare the results of Spectra, PageRank, and retweet count
indice=int64(Edge(:,1));
tally_edge = accumarray(indice, ones(size(indice)));
figure; hist(tally_edge, unique(tally_edge));
S_tally_edge=sort(tally_edge, 'descend');
topEdge=find(tally_edge>=S_tally_edge(10))
dlmwrite('top_edge.csv',[Index(topEdge) tally_edge(topEdge) Jt1(topEdge)], 'delimiter', ',', 'precision', 10);

% use <= for position dlmwrite('alternative_PR_result.csv', [PRScoreRank positionPR belowPR], 'delimiter', ',', 'precision', 10);

%% subgraph of ID
%% view number prediction
im=importdata('Count.txt');
Count=im.data;
names=im.textdata;

CT=diag(Count);
for i=1:length(Count)
    load(char(strcat(names(i), '.txt')));
end;
 A=[sum(BOAnews(:,2:4))
    ;sum(Citi(:,2:4))
    ;sum(MorganStanley(:,2:4))
    ;sum(WellsFargo(:,2:4))
    ;sum(Belk(:,2:4))
    ;sum(Lowes(:,2:4))
    ;sum(Macys(:,2:4))
    ;sum(Target(:,2:4))
    ;sum(Walmart(:,2:4))]
Am=[median(BOAnews(:,2:4))
    ;median(Citi(:,2:4))
    ;median(MorganStanley(:,2:4))
    ;median(WellsFargo(:,2:4))
    ;median(Belk(:,2:4))
    ;median(Lowes(:,2:4))
    ;median(Macys(:,2:4))
    ;median(Target(:,2:4))
    ;median(Walmart(:,2:4))];
CT*Am*0.71
CT*(A/150)

prediction=CT*Am*0.71;
dlmwrite('Prediction2.csv',prediction, 'delimiter', ',', 'precision', 10);
bank_retail_index=[9838356 1173193 4608767 142995 4198878 6120360 9635388 6856207 10886072];


%% Centrality
A1=A&A';
A2=double(logical(A1));
A2(~any(A2,2),:)=[];
A2(:,~any(A2,1))=[];
S=sum(A2);
S_d=S/(size(S,2)-1);
S1=1./S;
s=size(S,2);
p=(1:s);
W=sparse(p,p,S1,s,s);
WAW=W*A2*W+A2*W;
U=1+sum(WAW);
figure;
hist(U);
R=S';
L=logical(A2);
%Each state with path of length n+1
An=A2;
for i=1:5
    An=An*A2;
    LOA=logical(An);
    %get location of nodes connected with length n+1
    T=LOA-L;
    T=T>0;
    L=L|LOA;
    in=sum(double(T),2);
    R=[R in];
end;

sub=full(R(:,3));
s_b=sub(sub>0);
Z=accumarray(s_b, 1);
figure;
plot(Z);



%%Aug-adjcluster

