function [Sisters]=getSisters(clone)
%% gets Sister data for R
% uses Nordon clone format
Sisters.ProgenyID=[];
Sisters.Age=[];
Sisters.StopReason=[];
Sisters.Cause=[];
Sisters.DistanceToEdgeOfColony=[];
Sisters.ClusterID=[];
Sisters.Generation=[];
for i=1:length(clone)
    for j=1:length(clone{i}.track)
        if clone{i}.track{j}.TrackNum>1
            Sisters.ProgenyID(end+1)=clone{i}.track{j}.TrackNum;
            Sisters.Age(end+1)=clone{i}.track{j}.DeathTime-clone{i}.track{j}.BirthTime;
            Sisters.StopReason{length(Sisters.Age),1}=clone{i}.track{j}.StopReason;
            switch(clone{i}.track{j}.StopReason)
                case 'Division'
                    Sisters.Cause(end+1)=1;
                case 'Apoptosis'
                    Sisters.Cause(end+1)=2;
                case 'Endomitosis'
                    Sisters.Cause(end+1)=3;
                case 'Adhesion'
                    Sisters.Cause(end+1)=4;
                case 'Detachment'
                    Sisters.Cause(end+1)=5;
                case 'Not complete'
                    Sisters.Cause(end+1)=0;
                case 'Nil'
                    Sisters.Cause(end+1)=0;
                case 'Lost cell'
                    Sisters.Cause(end+1)=0;
                otherwise
                    Error('Fate not recognised');
            end
            Sisters.DistanceToEdgeOfColony(end+1)=mean(clone{i}.track{j}.DistanceToEdgeOfColony);
            Sisters.Generation(end+1)=floor(log2(clone{i}.track{j}.TrackNum));
            Sisters.ClusterID(end+1)=floor(clone{i}.track{j}.TrackNum/2)/128+i; % assign unique identify for sisters
        end   
    end
end
% assign cardinal instead of unique identifier, create lookup table 
ndx=1:length(Sisters.ClusterID);
ndx=[Sisters.ClusterID(:),ndx(:)];
ndx=sortrows(ndx,1);
ndx=[ndx,zeros(length(Sisters.ClusterID),1)]; % 3rd column is ClusterID
ndx(1,3)=1;
current=ndx(1,1);
for i=2:length(Sisters.ClusterID)
    if ndx(i,1)>current
        ndx(i,3)=ndx(i-1,3)+1;
        current=ndx(i,1);
    else
        ndx(i,3)=ndx(i-1,3);
    end
end
% 2nd column is ndx, 3rd column is ClusterID
for i=1:length(Sisters.ClusterID)
    Sisters.ClusterID(i)=ndx(i,3);
end
% convert to column vectors
Sisters.ProgenyID=Sisters.ProgenyID(:);
Sisters.Age=Sisters.Age(:);
Sisters.DistanceToEdgeOfColony=Sisters.DistanceToEdgeOfColony(:);
Sisters.Generation=Sisters.Generation(:);
Sisters.ClusterID=Sisters.ClusterID(:);
