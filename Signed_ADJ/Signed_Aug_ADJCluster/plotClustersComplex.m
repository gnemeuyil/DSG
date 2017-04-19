function [ output_args ] = plotClustersComplex( Jfs,d,vs ,I,N,M,Names )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% I: Column indicators for eigenpairs
% N: Row indicators for nodes
% M: Marker types for nodes
v=vs(:,I);
loc=(d(I)==real(d(I)));
if sum(~loc)>0
    v=[v(:,loc) real(v(:,~loc)) imag(v(:,~loc))];
end;
        
ps=idx2lgc(Jfs);
p=ps(:,I);
size(p)

p=logical(p);

pv=max(max(v))+0.1;
nv=min(min(v))-0.1;
dx = 0.03; dy = 0.03; dz=0.03;


Ind=false(size(Jfs,1),1);
for i=1:size(I,2)
    Ind=Ind|(Jfs==I(i));
end;

if size(v,2)>=3

figure;

hold on



plot3(v(p(:,1),1),v(p(:,1),2),v(p(:,1),3),'r^');

plot3(v(p(:,2),1),v(p(:,2),2),v(p(:,2),3),'b+');

%plot3(v(p(:,2),1),v(p(:,2),2),v(p(:,2),3),'go');



 
t=text(v(Ind,1)+dx, v(Ind,2)+dy,v(Ind,3)+dz, Names(Ind),'FontSize',16);
xlabel('R_1','FontSize',16);
ylabel('R_2','FontSize',16);
zlabel('I_2','FontSize',16);



axis([nv pv nv pv nv pv ]);







plot(0,0,'ko');
plot(0,0,'k.');
if ~isempty(N)
    for i=1:size(M,2)
        plot3(v(N(i),1),v(N(i),2),v(N(i),3),M{i});
    end
end;


%box on;
hold off;
daspect([1 1 1])


end



%if two dimentional
if size(v,2)<3
v=[v zeros(size(v,1),1)];
p=[p zeros(size(v,1),1)];
p=logical(p);



figure;

hold on



plot3(v(p(:,1),1),v(p(:,1),2),v(p(:,1),3),'r^');

plot3(v(p(:,2),1),v(p(:,2),2),v(p(:,2),3),'b+');

t=text(v(Ind,1)+dx, v(Ind,2)+dy,v(Ind,3)+dz, Names(Ind),'FontSize',16);


xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
%zlabel('x_3','FontSize',16);



axis([nv pv nv pv nv pv ]);







plot(0,0,'ko');
plot(0,0,'k.');
if ~isempty(N)
    for i=1:size(M,2)
        plot3(v(N(i),1),v(N(i),2),v(N(i),3),M{i});
    end
end;


%box on;
hold off;
daspect([1 1 1])



end


end

