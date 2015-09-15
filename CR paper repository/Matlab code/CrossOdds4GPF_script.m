%% Competition between absolute times of GFP expression and cell death
% Generates a data.frame for R with the following fields
%
% $Event_Time is the time of GFP onset taken from the beginning of the
% experiment
% 
% $Cause is the cause of failure, which in this case is 0=right censored,
% 1=GFP expression, 2=apoptosis
%
% $Clone is the clone ID
%
% $Growth_Factor is 1=G-CSF, 2=M-CSF
%
% $Progenitor is the progenitor ID
%
% $Generation is the generaiton

%% load g_clone and m_clone objects
load g_clone;
load m_clone;
g_lin=lineage(g_clone);
m_lin=lineage(m_clone);

%% Analysis by R functions in METS package: cor.cif, Grandom.cif, com.risk (timereg package)
header={'Event_Time' 'Start_Time' 'Cause' 'Growth_Factor' 'Clone' 'Progenitor' 'Generation'};
g_listmode=ShowData(g_lin,'Listmode');
m_listmode=ShowData(m_lin,'Listmode');
  % rows are ...
    % 1. clone number
    % 2. cell ID (track)
    % 3. cell generation
    % 4. birth phenotype         
    % 5. death phenotype
    % 6. number of phenotype transitions per lifetime
    % 7. transition age
    % 8. lifetime 
    % 9. stop reason None=0(right-censored), Division, Apoptosis, Lost
    %10. BirthAge

% find all branches with GFP transition
g_GFP_ndx=~isnan(g_listmode.wavelength(7,:));
m_GFP_ndx=~isnan(m_listmode.wavelength(7,:));
% look up GFP transition intervals
[g_GFP_time,g_GFP_StartTime]=getTransitionInterval(g_listmode.wavelength,g_GFP_ndx,'GFP_on');
[m_GFP_time,m_GFP_StartTime]=getTransitionInterval(m_listmode.wavelength,m_GFP_ndx,'GFP_on');
%% look up right censored intervals
g_RC_ndx=(g_listmode.wavelength(9,:)==0)&~...
    (getDescendants(g_listmode.wavelength,g_GFP_ndx,0)); % right censored
g_RC_ndx=g_RC_ndx|...
    ((g_listmode.wavelength(9,:)==3)&~...
    (getDescendants(g_listmode.wavelength,g_GFP_ndx,3))); % lost as well
m_RC_ndx=(m_listmode.wavelength(9,:)==0)&~...
    (getDescendants(m_listmode.wavelength,m_GFP_ndx,0)); % right censored
m_RC_ndx=m_RC_ndx|...
    ((m_listmode.wavelength(9,:)==3)&~...
    (getDescendants(m_listmode.wavelength,m_GFP_ndx,3))); % lost as well
[g_RC_time,g_RC_StartTime]=getTransitionInterval(g_listmode.wavelength,g_RC_ndx,'Right_Censored');
[m_RC_time,m_RC_StartTime]=getTransitionInterval(m_listmode.wavelength,m_RC_ndx,'Right_Censored');
%% look up apoptotic intervals
g_apop_ndx=(g_listmode.wavelength(9,:)==2)&~...
    (getDescendants(g_listmode.wavelength,g_GFP_ndx,2)); % apoptosis
m_apop_ndx=(m_listmode.wavelength(9,:)==2)&~...
    (getDescendants(m_listmode.wavelength,m_GFP_ndx,2)); % right censored
[g_apop_time,g_apop_StartTime]=getTransitionInterval(g_listmode.wavelength,g_apop_ndx,'Death');
[m_apop_time,m_apop_StartTime]=getTransitionInterval(m_listmode.wavelength,m_apop_ndx,'Death');


%% assemble fields
Event_Time=[g_GFP_time(:);m_GFP_time(:);g_RC_time(:);m_RC_time(:);g_apop_time(:);...
    m_apop_time(:)];
Start_Time=[g_GFP_StartTime(:);m_GFP_StartTime(:);g_RC_StartTime(:);m_RC_StartTime(:);...
    g_apop_StartTime(:);...
    m_apop_StartTime(:)];
Cause=[ones(size(g_GFP_time(:)));ones(size(m_GFP_time(:)));zeros(size(g_RC_time(:)));...
    zeros(size(m_RC_time(:)));2*ones(size(g_apop_time(:)));2*ones(size(m_apop_time(:)))];
Clone=[g_listmode.wavelength(1,g_GFP_ndx) m_listmode.wavelength(1,m_GFP_ndx) ...
    g_listmode.wavelength(1,g_RC_ndx) m_listmode.wavelength(1,m_RC_ndx)...
    g_listmode.wavelength(1,g_apop_ndx) m_listmode.wavelength(1,m_apop_ndx)];
Clone=Clone';
Progenitor=[g_listmode.wavelength(2,g_GFP_ndx) m_listmode.wavelength(2,m_GFP_ndx) ...
    g_listmode.wavelength(2,g_RC_ndx) m_listmode.wavelength(2,m_RC_ndx)...
    g_listmode.wavelength(2,g_apop_ndx) m_listmode.wavelength(2,m_apop_ndx)];
Progenitor=Progenitor';
Generation=[g_listmode.wavelength(3,g_GFP_ndx) m_listmode.wavelength(3,m_GFP_ndx) ...
    g_listmode.wavelength(3,g_RC_ndx) m_listmode.wavelength(3,m_RC_ndx)...
    g_listmode.wavelength(3,g_apop_ndx) m_listmode.wavelength(3,m_apop_ndx)];
Generation=Generation';
Growth_Factor=[ones(size(g_GFP_time(:)));2*ones(size(m_GFP_time(:)));ones(size(g_RC_time(:)));...
    2*ones(size(m_RC_time(:)));ones(size(g_apop_time(:)));2*ones(size(m_apop_time(:)))];

%% Save as text file for R
%'Event_Time' 'Cause' 'Growth_Factor' 'Clone' 'Progenitor' 'Generation'
% get rid of NaNs
b=~isnan(Event_Time);
Event_Time=Event_Time(b);
Cause=Cause(b);
Growth_Factor=Growth_Factor(b);
Clone=Clone(b);
Progenitor=Progenitor(b);
Generation=Generation(b);
Start_Time=Start_Time(b);
SavePrompt='Save competing hazards regression data for GFP onset';
Matlab2CSV4R(header,SavePrompt,Event_Time,Start_Time,Cause,Growth_Factor,Clone,Progenitor,Generation);

%% Next generate clonal pairs to analyse cross-odd between related cells
data=[Event_Time,Start_Time,Cause,Growth_Factor,Clone,Progenitor,Generation];
% sort with respect to Clone
data=sortrows(data,[4,5]);
% find clone IDs for G_CSF and M-CSF
G_ndx=data(:,4)==1;
M_ndx=data(:,4)==2;
G_Clone_IDs=unique(data(G_ndx,5));
M_Clone_IDs=unique(data(M_ndx,5));
% Generate clone pairs and assemble file for R
PairData=zeros(1,10);
PairID=1;
for i=1:length(G_Clone_IDs)
    b=data(:,5)==G_Clone_IDs(i);
    b=b&G_ndx;
    CloneData=data(b,:);
    [Pairs,Relatedness,Ancestor]=getPairs(CloneData(:,6)); % progeny pairs
    D=1e-10; % diffusion coefficient for small protein eg cytokine m2/s
    % uses the Green function (probability density) to calculate Likelihood
    % ratio for a 3 D diffusion process (infinite boundaries)
    % get T
    T=[];
    for j=1:length(Pairs(:))
        T(j)=CloneData((Pairs(j)==CloneData(:,6)),1)+CloneData((Pairs(j)==CloneData(:,6)),2);
    end
    if length(T)>1
        T=reshape(T,length(T)/2,2);
        [diff_t,del_x,LR]=getDiffusionTime(Pairs,T,G_Clone_IDs(i),D,g_clone,3,5); % 3-D cell radius 5 microns
    else
        diff_t=NaN;
        del_x=NaN;
        LR=NaN;
    end
    
    m=length(Relatedness);
    if m>0
        for j=1:m
            b1=Pairs(j,1)==CloneData(:,6);
            b2=Pairs(j,2)==CloneData(:,6);

            PairData=cat(1,PairData,[CloneData(b1,:) PairID Relatedness(j) LR(j);...
                CloneData(b2,:) PairID Relatedness(j) LR(j)]);
           
            PairID=PairID+1;
        end
    end
end
for i=1:length(M_Clone_IDs)
    b=data(:,5)==M_Clone_IDs(i);
    b=b&M_ndx;
    CloneData=data(b,:);
    [Pairs,Relatedness,Ancestor]=getPairs(CloneData(:,6)); % progeny pairs
    
    T=[];
    for j=1:length(Pairs(:))
        T(j)=CloneData((Pairs(j)==CloneData(:,6)),1)+CloneData((Pairs(j)==CloneData(:,6)),2);
    end
    if length(T)>1
        T=reshape(T,length(T)/2,2);
        [diff_t,del_x,LR]=getDiffusionTime(Pairs,T,M_Clone_IDs(i),D,m_clone,3,5); % 3-D cell radius 5 microns
    else
        diff_t=NaN;
        del_x=NaN;
        LR=NaN;
    end
    m=length(Relatedness);
    if m>0
        for j=1:m
            b1=Pairs(j,1)==CloneData(:,6);
            b2=Pairs(j,2)==CloneData(:,6);
                        PairData=cat(1,PairData,[CloneData(b1,:) PairID Relatedness(j) LR(j);...
                CloneData(b2,:) PairID Relatedness(j) LR(j)]);
            PairID=PairID+1;
        end
    end
end
header={'Event_Time' 'StartTime' 'Cause' 'Growth_Factor' 'Clone' 'Progenitor' 'Generation' 'PairID' 'Relatedness' 'Likelihood_Ratio'} ; 
PairData=PairData(2:end,:);
SavePrompt='Save Pair Data';
Matlab2CSV4R(header,SavePrompt,PairData(:,1),PairData(:,2),PairData(:,3),PairData(:,4),...
    PairData(:,5),PairData(:,6),PairData(:,7),PairData(:,8),PairData(:,9),PairData(:,10));



    




