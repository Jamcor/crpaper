function [Ancestor,Relationship,Relatedness]=getKinship(Progeny)
%% Calculates the Kinship code for two progeny
% Progeny is a 2 element vector with progeny numbers 1 is the ancestor, 2
% and 3 first daughters etc etc

% Ancestor is the common ancestor with the largest
% Relatedness is defined below
% Relationship is the text equivalent
% For example the Ancestor of progeny 2 and 3 is progeny 1. Ancestor 2 is
% the Ancestor for progeny 4 and 5, 8 and 5 etc.
% etc etc

%% Lookup table for text equivalent given (i,k)
% 	kth removed		abs(gen(Prog1)-gen(Prog2))			
% ith cousin	0	1	2	3	4	5
% 0	sisters	parents     grandparents	greatgrandparents	greatgreatgrandparents	
% 1	uncle	granduncle	greatgranduncle	greatgreatgranduncle	
% 2	1st cousin					
% 3	2nd cousin					
% 4	3rd cousin					
% 						
% 						
% I is min(Prog)-Gen(Anc)-1						
% 						
% 						
% Relatedness is k+1/(i+1)						
% 	k					
% i	0       1       2       3       4       5
% 0	1       2       3       4       5       6
% 1	0.5     1.5     2.5     3.5     4.5     5.5
% 2	0.3333  1.3333  2.3333  3.3333  4.3333  5.3333
% 3	0.25	1.25	2.25	3.25	4.25	5.25
% 4	0.2     1.2     2.2     3.2     4.2     5.2


%% find common ancestor
x=Progeny(1);
y=Progeny(2);
anc_x=x;
anc_y=y;
while x>1
    x=floor(x/2);
    anc_x=[x, anc_x];
end
while y>1
    y=floor(y/2);
    anc_y=[y, anc_y];
end
Ancestor=max(anc_x(ismember(anc_x,anc_y)));
GenAncestor=getGen(Ancestor);
GenProg=[getGen(Progeny(1)) getGen(Progeny(2))];
% are they sisters
if ((GenProg(1)-GenAncestor)==1)&&(GenProg(1)==GenProg(2))
    i=0;
    k=0;
% are they direct descendents?
elseif (Ancestor==Progeny(1))||(Ancestor==Progeny(2))
    i=0;
    k=abs(GenProg(1)-GenProg(2));
% are they uncles?
elseif (unique(min(GenProg))-GenAncestor)==1
    i=1;
    k=abs(GenProg(1)-GenProg(2));
% are they cousins
else
    k=abs(GenProg(1)-GenProg(2));
    i=unique(min(GenProg))-GenAncestor;
    %are they ancestors 

end
Relatedness=k+1/(i+1);   % same generation is 0-1, 1 generation difference 1-2...etc
Relationship=getRelationship(i,k);


    

    


function gen=getGen(x)
gen=floor(log2(x));

function Relationship=getRelationship(i,k)
% zeroth cousins
str={'sister' 'parent' 'grandparent'};
if (i==0)&&(k<3)
    Relationship=str{k+1};
elseif (i==0)&&(k>=3)
    Relationship=[repmat('great',1,k-2) 'grandparent'];
elseif (i==1)&&(k<3)
    Relationship=[repmat('grand',1,k-1) 'uncle'];
elseif (i==1)&&(k>=3)
    Relationship=[repmat('great',1,k-2) 'granduncle'];    
elseif (i>0)
    Relationship=[num2str(i-1) 'th cousin ' num2str(k) 'th removed'];
    
end


 
