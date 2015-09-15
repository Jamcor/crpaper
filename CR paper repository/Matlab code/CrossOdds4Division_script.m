%% Competition between absolute times of cell division and cell death
% Generates a data.frame for R with the following fields
%
% $Event_Time is the time lifespan of the cell 
% 
% $Cause is the cause of failure, which in this case is 0=right censored,
% 1=division, 2=apoptosis
%
% $Clone is the clone ID
%
% $Growth_Factor is 1=G-CSF, 2=M-CSF
%
% $GFP is the GFP status of the cell at death
%
% $Adherence is the adherence status of the cell (0, nonadherent, 1,
% adherent)
%
% $Progenitor is the progenitor ID
%
% $Generation is the generaiton

%% load g_clone and m_clone objects
load g_clone;
load m_clone;
g_lin=lineage(g_clone);
m_lin=lineage(m_clone);
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
    
% index generation>0
g_listmode=ShowData(g_lin,'Listmode');
m_listmode=ShowData(m_lin,'Listmode'); 
g_col=g_listmode.lineage(3,:)>0; % cell generation >0
m_col=m_listmode.lineage(3,:)>0;

% get Event_Time
g_Event_Time=g_listmode.lineage(8,g_col);
m_Event_Time=m_listmode.lineage(8,m_col);

% get GFP transition time
g_GFP_Age=g_listmode.lineage(7,g_col);
m_GFP_Age=m_listmode.lineage(7,m_col);

% get Cause
g_Cause=g_listmode.lineage(9,g_col);
m_Cause=m_listmode.lineage(9,m_col);

% get Growth Factor
g_GF=ones(size(g_Cause));
m_GF=2*ones(size(m_Cause));

% get GFP at end of lifetime
g_GFP=g_listmode.wavelength(5,g_col);
m_GFP=m_listmode.wavelength(5,m_col);

% get adherence at end of lifetime
g_adherent=g_listmode.lineage(5,g_col)>0;
m_adherent=m_listmode.lineage(5,m_col)>0;

% get generation
g_generation=g_listmode.lineage(3,g_col);
m_generation=m_listmode.lineage(3,m_col);

% get clone
g_cl=g_listmode.lineage(1,g_col);
m_cl=m_listmode.lineage(1,m_col);

% get Progenitor
g_progenitor=g_listmode.lineage(2,g_col);
m_progenitor=m_listmode.lineage(2,m_col);

% Assemble data.frame fields
Event_Time=[g_Event_Time';m_Event_Time'];
GFP_Age=[g_GFP_Age'; m_GFP_Age'];
Cause=[g_Cause';m_Cause'];
Growth_Factor=[g_GF';m_GF'];
Clone=[g_cl';m_cl'];
Progenitor=[g_progenitor';m_progenitor'];
Generation=[g_generation';m_generation'];
GFP=[g_GFP';m_GFP'];
Adherent=[g_adherent';m_adherent'];


%% Next generate clonal pairs to analyse cross-odd between related cells
data=[Event_Time,Cause,Growth_Factor,Clone,Progenitor,Generation,GFP,Adherent,GFP_Age];
% remove all records with NaN except GFP_Age
b=sum(isnan(data(:,1:8)),2)>0;
data=data(~b,:);
% find clone IDs for G_CSF and M-CSF
G_ndx=data(:,3)==1;
M_ndx=data(:,3)==2;
G_Clone_IDs=unique(data(G_ndx,4));
M_Clone_IDs=unique(data(M_ndx,4));
% Generate clone pairs and assemble file for R
PairData=zeros(1,14);
PairID=1;
for i=1:length(G_Clone_IDs)
    display(['Clone ' num2str(i) ' of ' num2str(length(G_Clone_IDs)) ' GCSF Clones']);
    b=data(:,4)==G_Clone_IDs(i);
    b=b&G_ndx;
    CloneData=data(b,:);
    [Pairs,Relatedness,Ancestor]=getPairs(CloneData(:,5)); % progeny pairs
    D=1e-10; % diffusion coefficient for small protein eg cytokine m2/s
    % uses the Green function (probability density) to calculate density in
    % 2 D after del_x and del_t with diffusion coefficient D for 2
    % dimensional infinite diffusion
    T=[];
    for j=1:length(Pairs(:))
        T(j)=CloneData((Pairs(j)==CloneData(:,5)),1);
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
            b1=Pairs(j,1)==CloneData(:,5);
            b2=Pairs(j,2)==CloneData(:,5);
            PairData=cat(1,PairData,[CloneData(b1,:) PairID Relatedness(j) diff_t(j) del_x(j) LR(j);...
                CloneData(b2,:) PairID Relatedness(j) diff_t(j) del_x(j) LR(j)]);
           
            PairID=PairID+1;
        end
    end
end
for i=1:length(M_Clone_IDs)
    display(['Clone ' num2str(i) ' of ' num2str(length(M_Clone_IDs)) ' MCSF Clones']);
    b=data(:,4)==M_Clone_IDs(i);
    b=b&M_ndx;
    CloneData=data(b,:);
    [Pairs,Relatedness,Ancestor]=getPairs(CloneData(:,5)); % progeny pairs
    T=[];
    for j=1:length(Pairs(:))
        T(j)=CloneData((Pairs(j)==CloneData(:,5)),1);
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
            b1=Pairs(j,1)==CloneData(:,5);
            b2=Pairs(j,2)==CloneData(:,5);
            PairData=cat(1,PairData,[CloneData(b1,:) PairID Relatedness(j) diff_t(j) del_x(j) LR(j);...
                CloneData(b2,:) PairID Relatedness(j) diff_t(j) del_x(j) LR(j)]);
            PairID=PairID+1;
        end
    end
end
header={'Event_Time' 'Cause' 'Growth_Factor' 'Clone' 'Progenitor' 'Generation' ...
    'GFP' 'Adherent' 'GFP_Age' 'PairID' 'Relatedness' 'Diffusion_time' 'Diffusion_distance' 'Likelihood_ratio'} ; 
PairData=PairData(2:end,:);
SavePrompt='Save Pair Data';
save('GMPPairData4Division','PairData','header');
Matlab2CSV4R(header,SavePrompt,PairData(:,1),PairData(:,2),PairData(:,3),PairData(:,4),...
    PairData(:,5),PairData(:,6),PairData(:,7),PairData(:,8),PairData(:,9),PairData(:,10),PairData(:,11),...
    PairData(:,12),PairData(:,13),PairData(:,14));
    




