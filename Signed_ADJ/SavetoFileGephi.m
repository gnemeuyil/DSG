function [ T_Nodes, T_Edges ] = SavetoFileGephi( subG_tilde,StartDate,S,NodeFile, EdgeFile)
%SAVETOFILEGEPHI Summary of this function goes here
%   Detailed explanation goes here
[Source,Target,Val]=find(subG_tilde);
P=Val>0;
N=~P;
Quote = repmat('''',nnz(Val),1);

% make the edge portion of gdf file
SDP=repmat(StartDate,sum(P),1);%positve start dates
SDN=(StartDate+1:StartDate+sum(N))';%negative start dates
SD=int2str([SDP;SDN]);%convert to strings
Quoted_SD=[Quote,SD,Quote];%append quotes
Weight=num2str(abs(Val),'%1.1f');
Quoted_Weight=[Quote,Weight,Quote];%convert weight to quoted strings as Edge width
EdgeLabel=num2str(sort(Val,'descend'),'%1.1f');% since Gephi don't take negative edge weights,
Quoted_EdgeLabel=[Quote,EdgeLabel,Quote];% use real edge weight as edge label
ColorB='''0,0,255''';
ColorR='''255,0,0''';
ColorBs=repmat(ColorB,sum(P),1);
ColorRs=repmat(ColorR,sum(N),1);
ColorEdges=[ColorBs;ColorRs];
SourceP=Source(P);
TargetP=Target(P);
[x,y]=ind2sub(size(subG_tilde),S);
Source=[SourceP;x'];
Target=[TargetP;y'];
Directed='true';
% make the edge table
T_Edges=table(Source,Target,repmat(Directed,nnz(Val),1),Quoted_EdgeLabel,Quoted_SD,Quoted_Weight,ColorEdges);

% make the node portion of gdf file
ColorG='''0,255,0''';
ColorNodes=repmat(ColorG,size(subG_tilde,1),1);
Node=1:size(subG_tilde,1);
QuoteNode = repmat('''',size(subG_tilde,1),1);
Node_Label=[QuoteNode,num2str(Node'),QuoteNode];
Vis=repmat('true',size(subG_tilde,1),1);
Label_Vis=repmat('true',size(subG_tilde,1),1);
% Make the node table
T_Nodes=table(Node',Node_Label,Vis,Label_Vis,ColorNodes);

% Write tables to two separate files that will be combined for Gephi input
writetable(T_Nodes,NodeFile,'Delimiter',',');
writetable(T_Edges,EdgeFile,'Delimiter',',');


end

