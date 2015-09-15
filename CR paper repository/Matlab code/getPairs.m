function [Pairs,Relatedness,Ancestor]=getPairs(Progeny)
%% Generates all possible relationships between Progeny which is a vector]
% Also calculates their relatedness
Pairs=combnk(Progeny,2);
n=size(Pairs,1);
Relatedness=zeros(n,1);
Ancestor=zeros(n,1);
for i=1:size(Pairs,1)
    [Ancestor(i),Relationship,Relatedness(i)]=getKinship(Pairs(i,:));
end