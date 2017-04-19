function [ G_7] = Combine7( G1,G2,G3,G4,G5,G6,G7,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
%COMBINE7 Summary of this function goes here
%   Detailed explanation goes here
% 0.1 0.2 0.1 0.1 0.15 0.15 0.2 0.1 0.1 0.1 0.1 0.1

%[ G_7,~,~] = CombineComponents( G1,G2,p1,p2);
[ G_7,~,~] = CombineComponents( G2,G3,p3,p4);
[ G_7,~,~] = CombineComponents( G_7,G4,p5,p6);
[ G_7,~,~] = CombineComponents( G_7,G5,p7,p8);
[ G_7,~,~] = CombineComponents( G_7,G6,p9,p10);
[ G_7,~,~] = CombineComponents( G_7,G7,p11,p12);
[ G_7,~,~] = CombineComponents( G1,G_7,p1,p2);
end

