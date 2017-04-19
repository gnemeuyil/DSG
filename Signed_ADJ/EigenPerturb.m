function [ D_subG_LR_tilde,V_subG_LR_tilde,D_subG_SR_tilde,V_subG_SR_tilde, subG_tilde ,S] = EigenPerturb( subG, K,val)
% Takes in a subgraph or cluster and perturb edges one by one.
% K is the 
[x,y]=find(subG);
I=sub2ind(size(subG),x,y);% convert edge subscript to index
Eye_ind=sub2ind(size(subG),1:size(subG,1),1:size(subG,1));% find the index of the diagonal elements
IND=1:size(subG,1)^2;% this is the index of all the possible edge entries
Edges = setdiff(IND,[I', Eye_ind]);% this is available edges indices without repeating the sub graph and diagonal elements


S = datasample(Edges,K,'Replace',false);%sample some without replacement to conduct the test

for i=1:K
    subG(S(i))=val;
    [v_subG_LR,d_subG_LR]=eigs(subG,ceil(size(subG,1)/2),'lr');
    DD_LR=diag(d_subG_LR);
    [~,I]=sort(real(DD_LR),'descend');
    D_subG_LR_tilde{i}=DD_LR(I);
    V_subG_LR_tilde{i}=v_subG_LR(:,I);
    %figure; hold on; plot(diag(d_subG_LR), 'r*');hold off
    [v_subG_SR,d_subG_SR]=eigs(subG,floor(size(subG,1)/2),'sr');
    DD_SR=diag(d_subG_SR);
    [~,I]=sort(real(DD_SR),'ascend');
    D_subG_SR_tilde{i}=DD_SR(I);
    V_subG_SR_tilde{i}=v_subG_SR(:,I);
end;
subG_tilde=subG;


end

