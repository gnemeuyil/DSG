function [me ma mi] = mdis(V,c,dist_type)
n = size(V,1);
C = ones(n,1)*c;


% Euclidean distance
if dist_type==1
    S = sqrt(sum((V - C).^2,2));
else
    % 1-cos distance
    S = 1-V*C';
end

me = mean(S);
ma = max(S);
mi = min(S);

