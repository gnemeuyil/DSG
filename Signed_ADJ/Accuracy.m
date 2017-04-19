function [ Acc ] = Accuracy( Jf0,Jf,md)
%ACCURACY Summary of this function goes here
%   Detailed explanation goes here
correct=0;
for i=1:max(Jf0{md})
%Idx=find(Jf0{md}==i);
correct=correct + sum(Jf{md}(Jf0{md}==i)==mode(Jf{md}(Jf0{md}==i)));
end
Acc=correct/size(Jf0{md},1);






end

