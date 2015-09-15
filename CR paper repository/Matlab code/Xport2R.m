function [ R_dataframe ] = Xport2R( clones, groups )
% Exports data from Matlab to R
%   clones is a cell array of clones
%   events are mutually exclusive competing risks e.g. mitosis, apoptos,
%   right censoring. String array
%   groups is structure array, first field is celltype, second field is treatment
%   and the second column are the treatments for each clone. The number of rows is equal to the numebr of clones
%   of clones and group are the same

    if length(clones)~=length(groups)
        error('Number of clones and number of rows in groups should be equal');
    end
    R_dataframe={'Age' 'Cell_Line' 'Treatment' 'StopReason' 'G1S_Age_1' 'G1S_Age_2' 'CloneID' 'TrackID' 'Generation'};
    for i=1:length(clones)
        R_dataframe=cat(1,R_dataframe,getCloneEvents(clones{i},groups(i)));
    end
   Matlab2CSV4R({R_dataframe{1,:}},'Save dataset',cell2mat({R_dataframe{2:end,1}}),...
       {R_dataframe{2:end,2}},...
       {R_dataframe{2:end,3}},...
       {R_dataframe{2:end,4}},...
       cell2mat({R_dataframe{2:end,5}}),...
       cell2mat({R_dataframe{2:end,6}}),...
       cell2mat({R_dataframe{2:end,7}}),...
       cell2mat({R_dataframe{2:end,8}}),...
       cell2mat({R_dataframe{2:end,9}}))
end

function Events=getCloneEvents(clone,group)
    LineNumber=1;
    for i=1:length(clone)
        for j=1:length(clone{i}.track)
            Events{LineNumber,1}=clone{i}.track{j}.DeathTime-clone{i}.track{j}.BirthTime; %Age
            Events{LineNumber,2}=group.celltype; %Cell_Line
            Events{LineNumber,3}=group.treatment; %Treatment
            Events{LineNumber,4}=clone{i}.track{j}.StopReason;  %StopReason          
            if length(clone{i}.track{j}.G1S)==1
                Events{LineNumber,5}=clone{i}.track{j}.G1S(1)-clone{i}.track{j}.BirthTime; %G1S_Age_1
                Events{LineNumber,6}=NaN;  %G1S_Age_2
            else
                Events{LineNumber,5}=clone{i}.track{j}.G1S(1)-clone{i}.track{j}.BirthTime;
                Events{LineNumber,6}=clone{i}.track{j}.G1S(2)-clone{i}.track{j}.BirthTime;
            end                
            Events{LineNumber,7}=i; %CloneID
            Events{LineNumber,8}=clone{i}.track{j}.TrackNum;%TrackID
            Events{LineNumber,9}=floor(log2(clone{i}.track{j}.TrackNum)); % generation
            LineNumber=LineNumber+1;
        end
    end
end




