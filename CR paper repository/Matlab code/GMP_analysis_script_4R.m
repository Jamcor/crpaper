%% Script for competing risk analysis in R
% Analysis of GMP data
% load g_clone and m_clone objects
load g_clone;
load m_clone;
g_lin=lineage(g_clone);
m_lin=lineage(m_clone);

%% CRA of mitosis and apoptosis age, without generation stratification, exclude generation 0
% fields for dataframe in R
% Event_Time Cause Growth_Factor GFP GCSF_GFPneg GCSF_GFPPos MCSF_GFPneg
% MCSF_GFPpos
header={'Event_Time' 'Cause' 'Growth_Factor' 'GFP' 'Group' 'Clone' 'Progeny' 'Generation' 'GFP_Age' 'IsGFPatBirth' 'StartTime' 'BirthTime'};
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

% index generation=>0
g_col=g_listmode.lineage(3,:)>=0;
m_col=m_listmode.lineage(3,:)>=0;

% get Event_Time
g_Event_Time=g_listmode.lineage(8,g_col);
m_Event_Time=m_listmode.lineage(8,m_col);
Event_Time=[g_Event_Time m_Event_Time];
% get Cause
g_Cause=g_listmode.lineage(9,g_col);
m_Cause=m_listmode.lineage(9,m_col);
Cause=[g_Cause m_Cause];
Cause(Cause==3)=0; % lost is also censored;
% get Clone
g_Clone=g_listmode.lineage(1,g_col);
m_Clone=m_listmode.lineage(1,m_col);
Clone=[g_Clone,m_Clone];
% get Progeny
g_Progeny=g_listmode.lineage(2,g_col);
m_Progeny=m_listmode.lineage(2,m_col);
Progeny=[g_Progeny,m_Progeny];
% get start of experiment and get cell birthtime
g_StartTime=zeros(size(g_Clone));
m_StartTime=zeros(size(m_Clone));
g_BirthTime=zeros(size(g_Clone));
m_BirthTime=zeros(size(m_Clone));
for i=1:length(g_StartTime)
    cloneID=g_Clone(i);
    progenyID=g_Progeny(i);
    ndx=(g_listmode.lineage(1,:)==cloneID)&(g_listmode.lineage(3,:)==0);
    g_StartTime(i)=g_listmode.lineage(10,ndx);
    ndx=(g_listmode.lineage(1,:)==cloneID)&(g_listmode.lineage(2,:)==progenyID);
    g_BirthTime(i)=g_listmode.lineage(10,ndx);
end
for i=1:length(m_StartTime)
    cloneID=m_Clone(i);
    progenyID=m_Progeny(i);
    ndx=(m_listmode.lineage(1,:)==cloneID)&(m_listmode.lineage(3,:)==0);
    m_StartTime(i)=m_listmode.lineage(10,ndx);
    ndx=(m_listmode.lineage(1,:)==cloneID)&(m_listmode.lineage(2,:)==progenyID);
    m_BirthTime(i)=m_listmode.lineage(10,ndx);
end
StartTime=[g_StartTime m_StartTime];
BirthTime=[g_BirthTime m_BirthTime];
    


% get Generation
g_Generation=g_listmode.lineage(3,g_col);
m_Generation=m_listmode.lineage(3,m_col);
Generation=[g_Generation,m_Generation];
% get Growth_Factor
Growth_Factor=cell(1,sum(g_col)+sum(m_col));
for i=1:sum(g_col)
    Growth_Factor{i}='gcsf';
end
for i=(sum(g_col)+1):(sum(g_col)+sum(m_col))
    Growth_Factor{i}='mcsf';
end
% get GFP (death phenotype)
g_GFP=g_listmode.wavelength(5,g_col);
m_GFP=m_listmode.wavelength(5,m_col);
GFP=[g_GFP m_GFP];
% isGFPatBirth
isg_GFPatBirth=g_listmode.wavelength(4,g_col);
ism_GFPatBirth=m_listmode.wavelength(4,m_col);
isGFPatBirth=[isg_GFPatBirth ism_GFPatBirth];
% get GPP age (transition age)
g_GFP_Age=g_listmode.wavelength(7,g_col);
m_GFP_Age=m_listmode.wavelength(7,m_col);
GFP_Age=[g_GFP_Age m_GFP_Age];

%get rid of columns with NaN
b=isnan(Event_Time)|isnan(Cause)|isnan(GFP);
Event_Time=Event_Time(~b);
Cause=Cause(~b);
Growth_Factor=Growth_Factor(~b);
GFP=GFP(~b);
Generation=Generation(~b);
Clone=Clone(~b);
GFP_Age=GFP_Age(~b);
isGFPatBirth=isGFPatBirth(~b);
StartTime=StartTime(~b);
BirthTime=BirthTime(~b);
Progeny=Progeny(~b)

% add 4 groups
Group=cell(size(Cause));
GCSF_GFPneg=(~GFP)&ismember(Growth_Factor,'gcsf');
GCSF_GFPpos=GFP&ismember(Growth_Factor,'gcsf');
MCSF_GFPneg=(~GFP)&ismember(Growth_Factor,'mcsf');
MCSF_GFPpos=GFP&ismember(Growth_Factor,'mcsf');
[Group{GCSF_GFPneg}]=deal('GCSF_GFPneg');
[Group{GCSF_GFPpos}]=deal('GCSF_GFPpos');
[Group{MCSF_GFPneg}]=deal('MCSF_GFPneg');
[Group{MCSF_GFPpos}]=deal('MCSF_GFPpos');
% create tab delimited file i.e. GMPregression
Matlab2CSV4R(header,'Save GMP regression data in R',Event_Time,Cause,Growth_Factor,GFP,...
    Group,Clone,Progeny,Generation,GFP_Age,isGFPatBirth,StartTime,BirthTime);
