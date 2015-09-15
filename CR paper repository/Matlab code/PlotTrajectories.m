function [dish] = PlotTrajectories( clones, listmode, outlier, max_time )
%% gets cell trajectories for all dishes using file field
% eg outlier = 1e4 microns
% max_time = 6 days
% plots positions of each clone as well as an ellipse including +/- 2 SD of
% the observations
% the dish structure include positions of GFP positive cells as well as the
% nearest neighbouring clone.
l=length(clones);
% get dish versus clone lookup table, unique dish ID is the date
dishlookup=zeros(l,3);
for i=1:l
    dishlookup(i,1)=i;
    str=clones{i}.header.FileName;
    start=strfind(str,'\');
    stop=strfind(str,'_');   
    dishlookup(i,2)=str2double(str((start(end)+1):stop(1)-1));
    if ~isempty(clones{i}.track{1})
        if ~isempty(clones{i}.track{1}.BirthTime)
            dishlookup(i,3)=clones{i}.track{1}.BirthTime;
        else
            dishlookup(i,3)=NaN;
        end
    end
end
dishlookup=sortrows(dishlookup,2);
dish_ids=unique(dishlookup(:,2));
map=jet(256);
for i=1:length(dish_ids)
    %get clone start time
    starttime=min(dishlookup(dish_ids(i)==dishlookup(:,2),3));
    %get clone IDs
    clone_ids=[];
    clone_ids=dishlookup(dishlookup(:,2)==dish_ids(i),1);
    
    dish(i).date=dish_ids(i);
    figure('name',['dish ' num2str(dish(i).date)]);
    for j=1:length(clone_ids)
        dish(i).clones{clone_ids(j),1}=clones{clone_ids(j)};
        GFPpositions=double(zeros(4,1)); % initialise GFP positions
        p=1;
        for k=1:length(clones{clone_ids(j)}.track)
            if ~isempty(clones{clone_ids(j)}.track{k})
                if ~isempty(clones{clone_ids(j)}.track{k}.X)                    
                    x=clones{clone_ids(j)}.track{k}.X(2,:);
                    y=clones{clone_ids(j)}.track{k}.Y(2,:);
                    t=clones{clone_ids(j)}.track{k}.Y(1,:);
                    b=(x<outlier)&(y<outlier);
                    x=x(b);
                    y=y(b);
                    t=t(b);

                    b=(listmode.wavelength(1,:)==clone_ids(j))&...
                        (listmode.wavelength(2,:)==clones{clone_ids(j)}.track{k}.TrackNum);                 

                   
                    if ~isnan(listmode.wavelength(7,b))
                        ndx=find(t>=(listmode.wavelength(7,b)+listmode.wavelength(10,b)),1,'first');
                        colour=map(ceil(256*(t(ndx)-starttime)/max_time),:);
                        plot(x(ndx),y(ndx),'.','MarkerEdgeColor',colour,'MarkerFaceColor',colour);
                        hold on;
%                         text(x(ndx),y(ndx),num2str(j),'FontSize',10);
                        %save cell position 
                        dish(i).clones{clone_ids(j)}.track{k}.GFPTransitionPosition=[x(ndx),y(ndx)];
                        dish(i).clones{clone_ids(j)}.track{k}.GFPTransitionTime=t(ndx);
                        if ~isempty(ndx)
                            GFPpositions(1,p)=clones{clone_ids(j)}.track{k}.TrackNum;
                            GFPpositions(2,p)=t(ndx);
                            GFPpositions(3,p)=x(ndx);
                            GFPpositions(4,p)=y(ndx); % first row is progenitor ID, second t, third x etc.                        
                            p=p+1;
                        end
                    end
                    dish(i).clones{clone_ids(j),1}.track{k}.X=[t;x];
                    dish(i).clones{clone_ids(j),1}.track{k}.Y=[t;y];                
                end
            end
        end       
         % plot centroid and 2SD region for each clone
        if p>1
            dish(i).clones{clone_ids(j)}.GFPTransitionPositions=GFPpositions;
            centroid=mean(GFPpositions(2:4,:),2);
            colour=map(ceil(256*(centroid(1)-starttime)/max_time),:);
            
%             plot(centroid(2),centroid(3),'*','MarkerSize',10,'MarkerEdgeColor',colour)
            sd=std(double(GFPpositions(3:4,:)),1,2);
            x=centroid(2)-2*sd(1);y=centroid(3)-2*sd(2); 
            w=4*sd(1);h=4*sd(2);
            if (w>0)&&(h>0)
                rectangle('position',[x,y,w,h],'curvature',[1,1],'EdgeColor',colour);
            end
            text(centroid(2),centroid(3),num2str(clone_ids(j)),'FontSize',10);
        else
            dish(i).clones{clone_ids(j)}.GFPTransitionPositions=[];
        end
         
    end
   
        
    % plot colour bar
    tickvalues=20:20:255;
    ticklabels=arrayfun(@(x) num2str(x,3),tickvalues*max_time/256,'UniformOutput',false);
    colormap(jet(256));
    colorbar('YTick',tickvalues,'YTickLabel',ticklabels);  
    dish(i).NearestClones=getNearestClone(dish(i));
end


function LookupTable=getNearestClone(dish)
Centroids=[];
for i=1:length(dish.clones)
    %get centroids
    if ~isempty(dish.clones{i})
        if ~isempty(dish.clones{i}.GFPTransitionPositions)
            v(1)=i;        
            v(2:3)=mean(dish.clones{i}.GFPTransitionPositions(3:4,:),2);
            Centroids=[Centroids; v];
        end
    end
end
n=size(Centroids,1);
LookupTable=zeros(n,2);
for i=1:n
    d2=sum((repmat(Centroids(i,2:3),[n,1])-Centroids(:,2:3)).^2,2);
    d2(i)=inf;
    b=d2==unique(min(d2));
    LookupTable(i,1)=Centroids(i,1);
    LookupTable(i,2)=Centroids(b,1);
end
% get rid of symmetric neighbours (counted twice)
z=LookupTable;
for i=1:n
    x=z(i,:);
    score=z-repmat([x(2) x(1)],[n,1]);
    score=sum((score(:,1)==0)&(score(:,2)==0));
    if score
        z(i,:)=NaN;
    end
    
end
LookupTable=LookupTable(~isnan(z(:,1)),:);
        
    
    

    
    
    