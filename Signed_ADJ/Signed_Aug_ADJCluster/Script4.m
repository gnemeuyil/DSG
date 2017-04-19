[ G_123,~,~] = CombineComponents( G3,G2,0,0);
[ G_123,~,~] = CombineComponents( G_123,G1,0,0);
denseity=sum(sum(G_123))/size(G_123,2)^2
[ Jf, Cf, kf, db, Q, acc] = Wcut( G_123,4,2,0,0);
[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v1_star,d1_star,v2_star,d2_star,I2,cand2] = Augmented_ADJCluster(G_123,4,2,0,0);
[Jf0, kf0, db0, Q0, ang0, acc0,Cf0,v1_0,d1_0,v2_0,d2_0,I0,cand0] = Augmented_ADJCluster(G_123,5,2,0,0);

[ G_123_02,~,~] = CombineComponents( G3,G2,0.2,0.2);
[ G_123_02,~,~] = CombineComponents( G_123_02,G1,0.2,0.2);


[ G_123_04,~,~] = CombineComponents( G3,G2,0.2,0.4);
[ G_123_04,~,~] = CombineComponents( G_123_04,G1,0.2,0.4);



[ G_123_04_05,~,~] = CombineComponents( G3,G2,0.4,0.5);
[ G_123_04_05,~,~] = CombineComponents( G_123_04_05,G1,0.4,0.5);


[ Jf, Cf, kf, db, Q, acc,V,D] = Wcut( G_123_04_05,10,2,0,0);
[Jf2, kf2, db2, Q2, ang2, acc2,Cf2,v1_star,d1_star,v2_star,d2_star,I2,cand2] = Augmented_ADJCluster(G_123_04_05,10,2,0,0);


%% Twitter4
loc=sum(A)<2;
sum(loc)

A1=A;
A1(loc,:)=[];
A1(:,loc)=[];
size(A1)


[Jfaa, kfaa, dbaa, Qaa, angaa, accaa,Cfaa,v1aa,d1aa,v2aa,d2aa,I2aa,cand2aa] = Augmented_ADJCluster(G_7_0_0_18,40,2,0,0);

dlmwrite('Jfaa.csv', Jfaa{2}, 'delimiter', ',', 'precision', 10);
dlmwrite('kfaa.csv', kfaa, 'delimiter', ',', 'precision', 10);
dlmwrite('dbaa.csv', dbaa, 'delimiter', ',', 'precision', 10);
dlmwrite('Qaa.csv', Qaa, 'delimiter', ',', 'precision', 10);
% the min mean and max angle of points to their centroids
dlmwrite('angaa.csv', angaa, 'delimiter', ',', 'precision', 10);
dlmwrite('accaa.csv', accaa, 'delimiter', ',', 'precision', 10);
dlmwrite('Cfaa.csv', Cfaa, 'delimiter', ',', 'precision', 10);
dlmwrite('v1aa.csv', v1aa, 'delimiter', ',', 'precision', 10);
dlmwrite('d1aa.csv', d1aa, 'delimiter', ',', 'precision', 10);
dlmwrite('v2aa.csv', v2aa, 'delimiter', ',', 'precision', 10);
dlmwrite('d2aa.csv', d2aa, 'delimiter', ',', 'precision', 10);
dlmwrite('I2aa.csv', I2aa, 'delimiter', ',', 'precision', 10);
dlmwrite('cand2aa.csv', cand2aa, 'delimiter', ',', 'precision', 10);


figure;plot(1:size(Qaa),abs(Qaa),'-o','MarkerSize',6,'MarkerFaceColor','b');
PI=[0,5, 6,9,12, 16, 17,18, 19,20, 21, 22, 24, 25];

figure;hold on;
plot(PI,abs(Qaa),'-','MarkerSize',6,'MarkerFaceColor','b');
plot(PI(1),abs(Qaa(1)),'^','MarkerSize',6,'MarkerFaceColor','b');
plot(PI([3, 6, 7, 9]),abs(Qaa([3, 6, 7, 9])),'o','MarkerSize',6,'MarkerFaceColor','r');
plot(PI([2, 4, 5, 8, 10, 11, 12, 13, 14]),abs(Qaa([2, 4, 5, 8, 10, 11, 12, 13, 14])),'s','MarkerSize',6,'MarkerFaceColor','g');
    set(gca,'XTick',PI(2:14))
   
    ylabel('Modularity','FontSize',14);
    hold off;
    box on;
   % axis square;


