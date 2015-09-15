function StartTimes=getAncestorStartTime(lin,clones)
% gets the start time of the experiment
% lin is the lineage in clone format (not an object)
% clones is a vector of clone IDs
StartTimes=zeros(size(clones));
for i=1:length(StartTimes)
    clone=clones(i);
    found=false;
    n=length(lin{clones(i)}.track);
    j=1;
    while ~found
        if lin{clone}.track{j}.TrackNum==1
            found=true;
            StartTimes(i)=lin{clone}.track{j}.BirthTime;
        else
            j=j+1;
            if j>n
                found=true;
                StartTimes(i)=NaN;
            end
        end
    end
end