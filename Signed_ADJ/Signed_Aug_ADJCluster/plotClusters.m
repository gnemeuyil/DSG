function [ output_args ] = plotClusters( Jfs,vs ,I,N,M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% I: Column indicators for eigenpairs
% N: Row indicators for nodes
% M: Marker types for nodes
v=vs(:,I);


ps=idx2lgc(Jfs);
p=ps(:,I);
size(p)

tic
p=logical(p);
td=0.01;
pv=max(max(v))+0.1;
nv=min(min(v))-0.1;

if size(v,2)>=3

figure;

hold on



plot3(v(p(:,1),1),v(p(:,1),2),v(p(:,1),3),'r^');

plot3(v(p(:,2),1),v(p(:,2),2),v(p(:,2),3),'b+');

plot3(v(p(:,3),1),v(p(:,3),2),v(p(:,3),3),'go');


xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
zlabel('x_3','FontSize',16);


%axis([nv pv nv pv nv pv ]);
axis([-0.6 0.6 -0.6 0.6 -0.6 0.6]);







plot(0,0,'ko');
plot(0,0,'k.');
if ~isempty(N)
    for i=1:size(M,2)
        plot3(v(N(i),1),v(N(i),2),v(N(i),3),M{i});
        t=text(v(N(i),1)+td,v(N(i),2)+td,v(N(i),3)+td, int2str(N(i)),'FontSize',16);
    end
end;


box on;
hold off;
daspect([1 1 1])
%box on
toc

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


xlabel('x_1','FontSize',16);
ylabel('x_2','FontSize',16);
zlabel('x_3','FontSize',16);


axis([-0.6 0.6 -0.6 0.6 -0.6 0.6]);
%axis([nv pv nv pv nv pv ]);







plot(0,0,'ko');
plot(0,0,'k.');
if ~isempty(N)
    for i=1:size(M,2)
        plot3(v(N(i),1),v(N(i),2),v(N(i),3),M{i});
        t=text(v(N(i),1)+td,v(N(i),2)+td,v(N(i),3)+td, int2str(N(i)),'FontSize',16);
    end
end;


box on;
hold off;
daspect([1 1 1])



end


end

