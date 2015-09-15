function ndx_descendents=getDescendants(listmode,ndx,fate)
%% gets terminal descendent of all cells index by ndx
ndx_descendents=false(size(ndx));
col=find(ndx);
for i=1:length(col)
    clone=listmode(1,col(i));
    maxgen=max(listmode(3,(listmode(1,:)==clone)));
    gen=listmode(3,col(i));
    proband=listmode(2,col(i));
    list=getlist(proband,gen,maxgen);% all descendents
    terminalbranch_ndx=listmode(9,:)==fate;% terminal branch events
    terminalbranch_ndx=(listmode(1,:)==clone)&terminalbranch_ndx; %terminal branches of clone    
    marked=ismember(list,listmode(2,terminalbranch_ndx));
    marked=list(marked);
    for j=1:length(marked)
        ndx_descendents=ndx_descendents|((listmode(2,:)==marked(j))&(listmode(1,:)==clone));
    end
end

function list=getlist(proband,gen,maxgen)
% get proband generation
list=proband;
nextgen=list;
for i=gen:(maxgen-1)    
    nextgen=[nextgen*2 nextgen*2+1]';
    list=[list;nextgen(:)];
end